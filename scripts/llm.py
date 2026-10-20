#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The language-model calls of the pipeline, isolated and cached, and
#           the readers that make them.
#           Backends (each reader names one; keys come from the environment only):
#             claude-cli   `claude -p` with every customization off (--safe-mode:
#                          no CLAUDE.md, skills, plugins, hooks, MCP servers), no
#                          tools, no session saved, auto-memory disabled, run
#                          from a scratch directory: the model sees the prompt
#                          alone. Empty model: the CLI's default model.
#             anthropic    the Anthropic Messages API through the official
#                          `anthropic` SDK (installed only for this backend);
#                          the key from the environment variable `api_key_env`
#                          (default ANTHROPIC_API_KEY) or the SDK's own profile.
#             openai       an OpenAI-compatible /chat/completions endpoint
#                          (`base_url`), the key from `api_key_env`.
#           Readers: two readers read the same records independently (neither sees
#           the other's answer); an adjudicator decides where they disagree. The
#           answer of every call must be one JSON value.
#           Cache: one file per (model, effort, system prompt, prompt) hash under
#           the cache directory (an answer cached under the CLI's default model also
#           serves the model that gave it, as the cache records it), holding the answer (e-mail addresses that are not
#           in the prompt removed), the model that answered and the time; a cached answer is reused. Records read in
#           batches are also cached one by one, so a record is read once whatever batch it falls in: a later snapshot
#           calls a model only for new or changed records. Usage and cost go to a local log
#           beside the cache (not kept in the repository). Offline (cache only):
#           no model is called; a question without a cached answer gets None.
#           Grounding helpers check an answer against the text it came from.
# ---------------------------------------------------------------------------
import hashlib
import http.client
import json
import os
import re
import subprocess
import threading
import time
import unicodedata
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

PROVIDERS = ('claude-cli', 'anthropic', 'openai')


class NoAnswer(Exception):
    """The model declined, or its answer was cut off: no answer for this question."""


class UsageLimit(RuntimeError):
    """The account's usage limit is reached: the run stops. Every answer received so far is cached, so a run after
    the limit resets continues where this one stopped."""


LIMIT_WORDS = ('usage limit', 'rate limit', 'limit reached', 'quota', 'credit balance', 'too many requests')


class Reader:
    """One model behind one backend, with the shared answer cache."""

    def __init__(self, cache_dir, workdir, model='', effort='', timeout=600, parallel=3, provider='claude-cli',
                 bin='claude', api_key_env='', base_url='', max_tokens=16000, offline=False, name='A'):
        assert provider in PROVIDERS, f'unknown provider {provider}; one of {PROVIDERS}'
        assert provider == 'claude-cli' or model, f'reader {name}: the {provider} backend needs a model'
        self.cache = Path(cache_dir)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.name, self.provider, self.model, self.effort = name, provider, model or '', effort or ''
        self.timeout, self.parallel, self.bin, self.max_tokens = timeout, parallel, bin or 'claude', max_tokens
        self.api_key_env, self.base_url, self.offline = api_key_env, base_url.rstrip('/'), offline
        self.calls, self.cached, self.missed, self.failed = 0, 0, 0, 0
        self.models_seen = set()
        self._lock = threading.Lock()
        self._client = None

    # the cache key: unchanged for the CLI and the Anthropic API (the answers of a model are its answers); an
    # OpenAI-compatible model is keyed with its endpoint, so that names of different services never collide
    def key(self, system, prompt):
        mk = self.model if self.provider in ('claude-cli', 'anthropic') else f'openai:{self.base_url}:{self.model}'
        h = hashlib.sha256('\x1e'.join([mk, self.effort, system, prompt]).encode()).hexdigest()
        return self.cache / h[:2] / f'{h}.json'

    def record_key(self, system, record):
        """The cache file of one record's answer, read in a batch (per_record)."""
        return self.key(system, 'One record of a batch:\n\n' + record_text(record))

    def record_cached(self, system, record):
        """The cached answer of one record read in an earlier batch, or None."""
        f = self.record_key(system, record)
        if not f.exists():
            return None
        rec = json.loads(f.read_text())
        with self._lock:
            self.cached += 1
            self.models_seen.update(rec.get('model') or [])
        return rec['answer']

    def cache_record(self, system, record, answer, models):
        """Cache one record's answer, taken from the answer of its batch."""
        f = self.record_key(system, record)
        if not f.exists():
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_text(json.dumps({'answer': answer, 'model': models, 'at': time.strftime('%F %T'),
                                     'from': 'batch'}, ensure_ascii=False))

    def describe(self):
        return f'{self.provider} {self.model or "(default model)"}' + (f', effort {self.effort}' if self.effort else '')

    def default_answer(self, system, prompt):
        """A cached answer of the command-line backend's default model (cached under an empty model name) that was
        given by this reader's model, as the cache records it: the answers of the first CTEPH build were cached so,
        by claude-opus-5. None if there is none."""
        if self.provider != 'claude-cli' or not self.model:
            return None
        h = hashlib.sha256('\x1e'.join(['', self.effort, system, prompt]).encode()).hexdigest()
        f = self.cache / h[:2] / f'{h}.json'
        if not f.exists():
            return None
        rec = json.loads(f.read_text())
        return rec if any(re.sub(r'\[.*\]$', '', m) == self.model for m in rec.get('model') or []) else None

    def ask(self, system, prompt, tries=3):
        """One answer (a parsed JSON value), from the cache or the model; None offline without a cached answer,
        when the model declines, or when its answer is no JSON after `tries` calls."""
        return self.ask_rec(system, prompt, tries)[0]

    def ask_rec(self, system, prompt, tries=3):
        """(answer, the models that gave it); see ask. A backend that keeps failing (network, CLI) stops the run:
        every answer received so far is cached."""
        f = self.key(system, prompt)
        rec = json.loads(f.read_text()) if f.exists() else self.default_answer(system, prompt)
        if rec is not None:
            with self._lock:
                self.cached += 1
                self.models_seen.update(rec.get('model') or [])
            return rec['answer'], rec.get('model') or []
        if self.offline:
            with self._lock:
                self.missed += 1
            return None, []
        err, malformed = '', False
        for k in range(tries):
            try:
                text, models, usage = self._call(system, prompt)
            except NoAnswer as e:
                with self._lock:
                    self.failed += 1
                self._usage({'reader': self.name, 'event': 'no answer', 'detail': str(e)[:200]})
                return None, []
            except UsageLimit:
                raise
            except (KeyError, RuntimeError) as e:
                err, malformed = str(e)[:400], False
                time.sleep(5 * (k + 1))
                continue
            try:
                answer = parse_json(text)
            except ValueError as e:                 # the model answered, but not with JSON
                err, malformed = str(e)[:400], True
                time.sleep(5 * (k + 1))
                continue
            answer, removed = scrub(answer, system + prompt)
            if removed:
                self._usage({'reader': self.name, 'event': 'e-mail addresses removed from an answer', 'n': removed})
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_text(json.dumps({'answer': answer, 'model': models, 'at': time.strftime('%F %T')},
                                    ensure_ascii=False))
            with self._lock:
                self.calls += 1
                self.models_seen.update(models)
            self._usage({'reader': self.name, 'event': 'call', 'model': models, **usage})
            return answer, models
        if malformed:                               # this question gets no answer; the run goes on
            with self._lock:
                self.failed += 1
            self._usage({'reader': self.name, 'event': 'no answer', 'detail': f'no JSON after {tries} calls: {err}'})
            return None, []
        raise RuntimeError(f'reader {self.name}: no answer after {tries} tries: {err}')

    def map(self, jobs, full=False):
        """[(system, prompt), ...] -> answers in order (full: (answer, models)), `parallel` calls at a time."""
        with ThreadPoolExecutor(self.parallel) as ex:
            return list(ex.map(lambda j: self.ask_rec(*j) if full else self.ask(*j), jobs))

    def stats(self):
        return {'reader': self.name, 'backend': self.describe(), 'calls': self.calls, 'cached': self.cached,
                'not_generated': self.missed, 'no_answer': self.failed, 'models': sorted(self.models_seen)}

    # ── backends ─────────────────────────────────────────────────────────────
    def _call(self, system, prompt):
        if self.provider == 'claude-cli':
            return self._claude_cli(system, prompt)
        if self.provider == 'anthropic':
            return self._anthropic(system, prompt)
        return self._openai(system, prompt)

    def _claude_cli(self, system, prompt):
        cmd = [self.bin, '-p', '--safe-mode', '--tools', '', '--strict-mcp-config', '--no-session-persistence',
               '--output-format', 'json', '--system-prompt', system]
        if self.model:
            cmd += ['--model', self.model]
        if self.effort:
            cmd += ['--effort', self.effort]
        env = {**os.environ, 'CLAUDE_CODE_DISABLE_AUTO_MEMORY': '1'}
        env.pop('CLAUDECODE', None)
        r = subprocess.run(cmd, input=prompt, capture_output=True, text=True, cwd=self.workdir, env=env,
                           timeout=self.timeout)
        try:
            out = json.loads(r.stdout)
        except json.JSONDecodeError:
            raise RuntimeError(f'no JSON from the CLI; stderr: {r.stderr[-300:]}')
        if out.get('is_error'):
            msg = str(out.get('result', ''))[:300]
            if any(w in msg.lower() for w in LIMIT_WORDS):
                raise UsageLimit(f'reader {self.name}: {msg}')
            raise RuntimeError(msg)
        usage = {'cost_usd': out.get('total_cost_usd'), 'duration_ms': out.get('duration_ms')}
        return out['result'], sorted((out.get('modelUsage') or {}).keys()), usage

    def _anthropic(self, system, prompt):
        import anthropic           # only this backend needs the SDK
        if self._client is None:
            key = os.environ.get(self.api_key_env or 'ANTHROPIC_API_KEY')
            self._client = anthropic.Anthropic(api_key=key, max_retries=5, timeout=self.timeout) if key else \
                anthropic.Anthropic(max_retries=5, timeout=self.timeout)
        kw = {'model': self.model, 'max_tokens': self.max_tokens, 'system': system,
              'messages': [{'role': 'user', 'content': prompt}]}
        if self.effort:
            kw['output_config'] = {'effort': self.effort}
        try:
            r = self._client.messages.create(**kw)
        except anthropic.RateLimitError as e:
            raise UsageLimit(f'reader {self.name}: rate limit after the SDK\'s retries: {e}')
        except (anthropic.InternalServerError, anthropic.APIConnectionError) as e:
            raise RuntimeError(f'{type(e).__name__}: {e}')
        if r.stop_reason == 'refusal':
            raise NoAnswer(f'declined ({getattr(r.stop_details, "category", None)})')
        if r.stop_reason == 'max_tokens':
            raise NoAnswer('answer cut off at max_tokens')
        text = ''.join(b.text for b in r.content if b.type == 'text')
        u = r.usage
        return text, [r.model], {'input_tokens': u.input_tokens, 'output_tokens': u.output_tokens}

    def _openai(self, system, prompt):
        key = os.environ.get(self.api_key_env or 'OPENAI_API_KEY', '')
        body = json.dumps({'model': self.model, 'messages': [{'role': 'system', 'content': system},
                                                             {'role': 'user', 'content': prompt}]}).encode()
        headers = {'Content-Type': 'application/json', **({'Authorization': f'Bearer {key}'} if key else {})}
        req = urllib.request.Request(f'{self.base_url}/chat/completions', data=body, headers=headers)
        try:
            d = json.load(urllib.request.urlopen(req, timeout=self.timeout))
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504, 529):
                raise RuntimeError(f'HTTP {e.code}')
            raise NoAnswer(f'HTTP {e.code}: {e.read().decode()[:200]}')
        except (urllib.error.URLError, TimeoutError, ConnectionError, http.client.HTTPException,
                json.JSONDecodeError) as e:                    # a network failure, or a response cut short: retried
            raise RuntimeError(str(e))
        ch = (d.get('choices') or [{}])[0]
        if ch.get('finish_reason') == 'length':
            raise NoAnswer('answer cut off at the length limit')
        u = d.get('usage') or {}
        return (ch.get('message') or {}).get('content') or '', [d.get('model') or self.model], \
            {'input_tokens': u.get('prompt_tokens'), 'output_tokens': u.get('completion_tokens')}

    def _usage(self, rec):
        """Usage and cost of the calls, in a local log beside the cache (git-ignored)."""
        try:
            with self._lock, open(self.cache.parent / 'llm_usage.jsonl', 'a') as fh:
                fh.write(json.dumps({'at': time.strftime('%F %T'), **rec}, ensure_ascii=False) + '\n')
        except OSError:
            pass


def load_llm(path, offline=False):
    """(readers, adjudicator) from the run's model configuration (llm.json, written by main.nf from params.llm):
    {readers: [{name, provider, model, effort, bin, api_key_env, base_url}, ...], adjudicator: {...} | null,
     cache_dir, workdir, parallel, timeout, max_tokens}. No key is ever in this file, only the names of the
    environment variables that hold them."""
    c = json.loads(Path(path).read_text())
    common = {'cache_dir': c['cache_dir'], 'workdir': c['workdir'], 'parallel': int(c.get('parallel') or 3),
              'timeout': int(c.get('timeout') or 600), 'max_tokens': int(c.get('max_tokens') or 16000),
              'offline': offline}
    fields = ('provider', 'model', 'effort', 'bin', 'api_key_env', 'base_url')

    def make(spec, name):
        return Reader(name=name, **common, **{k: spec[k] for k in fields if spec.get(k) is not None})
    readers = [make(r, r.get('name') or 'AB'[i]) for i, r in enumerate(c['readers'])]
    assert 1 <= len(readers) <= 2, 'one or two readers'
    adj = c.get('adjudicator')
    assert len(readers) == 1 or adj, 'two readers need an adjudicator: without one, every record the readers ' \
                                     'disagree on would be left out'
    # the adjudicator: its own model and effort; the backend of reader A unless it names one
    backend = {k: c['readers'][0][k] for k in ('provider', 'bin', 'api_key_env', 'base_url') if k in c['readers'][0]}
    adjudicator = make({**backend, **adj}, 'adjudicator') if adj and len(readers) == 2 else None
    if adjudicator:
        same = [r.name for r in readers if (r.provider, r.model) == (adjudicator.provider, adjudicator.model)]
        assert not same, f'the adjudicator ({adjudicator.describe()}) must be a third model, not that of reader {same}'
    return readers, adjudicator


def read_all(readers, jobs):
    """Every reader answers every job, the readers side by side: [[answers of reader A], [answers of reader B]]."""
    if len(readers) == 1:
        return [readers[0].map(jobs)]
    with ThreadPoolExecutor(len(readers)) as ex:
        return list(ex.map(lambda r: r.map(jobs), readers))


EMAIL = re.compile(r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}')


def scrub(answer, given):
    """(answer, n): the answer without the e-mail addresses that are not in the text it was given. A command-line
    backend logged in to an account may add the account's address to what the model sees; it must never reach an
    answer, which is cached and may be published."""
    n = [0]

    def walk(x):
        if isinstance(x, str):
            def rep(m):
                if m.group(0) in given:
                    return m.group(0)
                n[0] += 1
                return '<e-mail address removed>'
            return EMAIL.sub(rep, x)
        if isinstance(x, list):
            return [walk(v) for v in x]
        if isinstance(x, dict):
            return {k: walk(v) for k, v in x.items()}
        return x
    return walk(answer), n[0]


def parse_json(text):
    """The JSON value of an answer, with or without a ```json fence. Text after the JSON is ignored, and so is
    broken JSON after a complete value. Several values in a row (with or without commas between them): arrays are
    joined; objects that each name a record (`pmid`) are the list of them; other objects are merged into one (lists
    joined, else the first value), so that no item a model listed in a second object is lost."""
    t = text.strip()
    m = re.search(r'```(?:json)?\s*(.*?)```', t, re.S)
    if m:
        t = m.group(1).strip()
    if not t.startswith(('{', '[')):
        i = min([x for x in (t.find('{'), t.find('[')) if x >= 0], default=-1)
        if i < 0:
            raise ValueError('no JSON in the answer')
        t = t[i:]
    dec, values, i = json.JSONDecoder(), [], 0
    while i < len(t) and t[i] in '{[':
        try:
            v, i = dec.raw_decode(t, i)
        except json.JSONDecodeError:
            if values:                      # broken text after a complete value
                break
            raise
        values.append(v)
        while i < len(t) and (t[i].isspace() or t[i] == ','):
            i += 1
    if len(values) == 1:
        return values[0]
    if all(isinstance(v, list) for v in values):
        return [x for v in values for x in v]
    if not all(isinstance(v, dict) for v in values):
        return values[0]
    if all('pmid' in v for v in values):
        return values
    merged = {}
    for v in values:
        for k, x in v.items():
            if k not in merged:
                merged[k] = list(x) if isinstance(x, list) else x
            elif isinstance(merged[k], list) and isinstance(x, list):
                merged[k] += x
    return merged


# ── grounding ────────────────────────────────────────────────────────────────
def unify(s):
    """Text for comparison: Unicode-normalised, dashes, quotes and spaces unified (case kept)."""
    s = unicodedata.normalize('NFKC', str(s or ''))
    s = re.sub(r'[‐-―−]', '-', s)
    s = s.replace(' ', ' ').replace('’', "'").replace('“', '"').replace('”', '"')
    return ' '.join(s.split())


def norm(s):
    """Text for substring checks: unify(), lower case."""
    return unify(s).lower()


def in_words(word, text, case=True):
    """word occurs in text as a whole word (not inside a longer word or number), in its own case unless case=False:
    'Asian' is not in 'Caucasian', 'China' not in 'Chinatown', 'F5' not in 'F50'."""
    w, t = unify(word), unify(text)
    return bool(w) and re.search(rf'(?<![A-Za-z0-9]){re.escape(w)}(?![A-Za-z0-9])', t, 0 if case else re.I) is not None


# places that put a study in East Asia when a text names them (with the country and population words of
# ot_common.EAS_WORDS): a hospital "in Peking" places the study as surely as "Chinese patients" does
EAST_ASIA_PLACES = ('Beijing', 'Peking', 'Shanghai', 'Guangzhou', 'Wuhan', 'Chengdu', 'Nanjing', 'Hangzhou',
                    'Taipei', 'Tokyo', 'Osaka', 'Nagoya', 'Seoul', 'Busan')


def in_symbol(name, text):
    """A gene or variant name occurs in text: a short one (at most 4 characters, e.g. F5, AR, ACE) as a whole word in
    its own case, a longer one as in_text (a short name inside other words is no evidence: 'one of 5 patients')."""
    return in_words(name, text) if len(re.sub(r'\s', '', str(name or ''))) <= 4 else in_text(name, text)


def in_text(fragment, text):
    """fragment occurs in text (normalised; for gene / variant names also without '-', '*' and ':')."""
    f, t = norm(fragment), norm(text)
    if not f:
        return False
    if f in t:
        return True
    strip = lambda x: re.sub(r'[-*:()\[\]\s]', '', x)  # noqa: E731   PubMed writes HLA alleles as DPB1(*)0202
    return strip(f) in strip(t)


def numbers_in(value, text):
    """Every number of value (e.g. '1,901' or 1901) occurs in text."""
    t = norm(text).replace(',', '')
    for n in re.findall(r'\d+(?:\.\d+)?', str(value).replace(',', '')):
        if not re.search(rf'(?<![\d.]){re.escape(n)}(?![\d])', t):
            return False
    return True


def overlap(a, b):
    """Two quotes of one text agree: one contains the other, or they share most of their words."""
    a, b = norm(a), norm(b)
    if not a or not b:
        return False
    if a in b or b in a:
        return True
    wa, wb = set(re.findall(r'\w+', a)), set(re.findall(r'\w+', b))
    return len(wa & wb) >= 0.6 * min(len(wa), len(wb))


# ── two readers, one adjudicator ─────────────────────────────────────────────
def record_text(r):
    """A PubMed record as the readers are given it."""
    return (f'PMID: {r["pmid"]}\nTitle: {r["title"]}\nPublication types: {"; ".join(r["pub_types"])}\n'
            f'Abstract: {r["abstract"] or "(no abstract)"}')


ADJUDICATE = ('\n\nYou are given one record and two readings of it, made by two readers who did not see each other\'s '
              'answer. Where the readings differ, decide from the record text alone which is right, or that neither '
              'is, and return your own reading of the record in the format above, with one more field: "decision": '
              '"<at most 25 words: what you decided where the readings differed, and why>".')


def per_record(readers, recs, system, batch, ground, logs):
    """Every reader reads the records, `batch` per call: {pmid: [grounded reading of each reader, or None]}.
    ground(answer of one record, record, reader name) -> the grounded reading. A record whose own answer is cached
    (from an earlier batch) is not read again; the others are read in batches, in their order, and each answer is
    cached for its record."""
    def read(rd):
        own = {r['pmid']: rd.record_cached(system, r) for r in recs}
        todo = [r for r in recs if own[r['pmid']] is None]
        bs = [todo[i:i + batch] for i in range(0, len(todo), batch)]
        answers = rd.map([(system, 'Records:\n\n' + '\n\n---\n\n'.join(record_text(r) for r in b)) for b in bs],
                         full=True)
        for b, (ans, models) in zip(bs, answers):
            ans = ans if isinstance(ans, list) else [ans] if isinstance(ans, dict) else []
            by = {str(x.get('pmid')): x for x in ans if isinstance(x, dict)}
            for r in b:
                if r['pmid'] in by:
                    own[r['pmid']] = by[r['pmid']]
                    rd.cache_record(system, r, by[r['pmid']], models)
        return own
    with ThreadPoolExecutor(len(readers)) as ex:
        answers = list(ex.map(read, readers))
    readings = {r['pmid']: [None] * len(readers) for r in recs}
    for ri, rd in enumerate(readers):
        for r in recs:
            x = answers[ri][r['pmid']]
            if isinstance(x, dict):
                readings[r['pmid']][ri] = ground(x, r, rd.name)
    for r in recs:
        if all(x is None for x in readings[r['pmid']]):
            logs.append({'pmid': r['pmid'], 'item': 'record', 'value': 'no answer', 'reader': ''})
    return readings


def reconcile(readers, adjudicator, items, system, differences, ground):
    """The final reading of each item read by the readers. items: [(key, source text, [readings])];
    differences(reading A, reading B) -> [(item, value A, value B)] of the decisive items on which they differ;
    ground(answer, key, source text) -> the grounded reading of the adjudicator's answer.
    One reader: its reading (agreement 'one reader'). Two readers that agree: reader A's reading ('both'). Otherwise
    the adjudicator reads the source text with both readings ('adjudicated'); a reader without an answer counts as a
    difference. Returns ({key: (final reading, agreement)}, review rows [{key, item, A, B, decision}])."""
    finals, review, jobs, todo = {}, [], [], []
    names = [r.name for r in readers]
    for key, text, rs in items:
        if all(x is None for x in rs):
            continue
        if len(readers) == 1:
            finals[key] = (rs[0], 'one reader')
            continue
        diff = differences(*rs) if None not in rs else [('reading', 'answered' if rs[0] is not None else 'no answer',
                                                         'answered' if rs[1] is not None else 'no answer')]
        if not diff:
            finals[key] = (rs[0], 'both')
            continue
        review += [{'key': key, 'item': it, names[0]: va, names[1]: vb} for it, va, vb in diff]
        # the readings in an order that alternates by record (fixed by its key), so that neither reader is
        # always the first the adjudicator reads
        swapped = bool(int(hashlib.sha256(str(key).encode()).hexdigest(), 16) % 2)
        first, second = (rs[1], rs[0]) if swapped else rs
        jobs.append((system + ADJUDICATE, f'{text}\n\n---\n\nReading 1:\n'
                     + json.dumps(first if first is not None else 'no answer', ensure_ascii=False)
                     + '\n\nReading 2:\n' + json.dumps(second if second is not None else 'no answer',
                                                       ensure_ascii=False)))
        todo.append((key, text, swapped))
    answers = adjudicator.map(jobs) if adjudicator and jobs else [None] * len(jobs)
    for (key, text, swapped), ans in zip(todo, answers):
        ans = ans[0] if isinstance(ans, list) and ans else ans
        if not isinstance(ans, dict):
            for row in review:
                if row['key'] == key:
                    row['decision'] = 'no answer from the adjudicator: left out'
            continue
        finals[key] = (ground(ans, key, text), 'adjudicated')
        for row in review:
            if row['key'] == key:
                row['decision'] = reader_names(str(ans.get('decision') or '')[:300], swapped, names)
    return finals, review


def decide(readers, adjudicator, system, batches, ground, show, note=''):
    """Decisions on a list of items, made the way every reading is: each reader decides every item without seeing the
    other's answer, the decisions are compared item by item, and the adjudicator decides the items where they differ.
    batches: [(text, ids)], the text the readers are given for the items `ids` (a long list is read in parts);
    ground(answer, ids, who) -> {id: decision}, the valid decisions of an answer, each a tuple whose last element is
    a reason (not compared); show(decision) -> its text, for the adjudicator and the review table; note: what the
    adjudicator is asked beyond the readers' instructions. An item no reader decided is left out of the result (the
    caller keeps it as it was). Returns ({id: decision}, {id: agreement}, review rows [{item, <reader>..., final,
    decision}], answered: whether any reader answered)."""
    names = [r.name for r in readers]
    answers = read_all(readers, [(system, text) for text, _ids in batches])
    if all(a is None for rd in answers for a in rd):
        return {}, {}, [], False
    final, agree, review = {}, {}, []
    for b, (text, ids) in enumerate(batches):
        reads = [ground(answers[k][b], ids, names[k]) for k in range(len(readers))]
        disputed = []
        for i in ids:
            ds = [rd.get(i) for rd in reads]
            if all(d is None for d in ds):
                continue
            if len(ds) == 1 or (None not in ds and ds[0][:-1] == ds[1][:-1]):
                final[i], agree[i] = ds[0], 'both' if len(ds) > 1 else 'one reader'
            else:
                disputed.append(i)
        if not disputed:
            continue
        # the readings in an order fixed by the part's text, so that neither reader is always shown first
        swapped = bool(int(hashlib.sha256(text.encode()).hexdigest(), 16) % 2)
        lines = []
        for i in disputed:
            a, c = (reads[k].get(i) for k in ((1, 0) if swapped else (0, 1)))
            lines.append(f'[{i}] Reading 1: {show(a) if a else "no decision"}; Reading 2: {show(c) if c else "no decision"}')
        ans = adjudicator.ask(system, text + '\n\n---\n\nTwo readers, who did not see each other\'s answer, decided '
                              'differently on these items:\n' + '\n'.join(lines) + '\n\nDecide each of these items from '
                              'the text above, in the format above, one decision per item listed here; each reason says '
                              'in at most 25 words what you decided where the readings differed, and why.'
                              + (f' {note}' if note else '')) if adjudicator else None
        ad = ground(ans, disputed, 'adjudicator')
        for i in disputed:
            if i in ad:
                final[i] = ad[i]
            agree[i] = 'adjudicated'
            review.append({'item': i, **{n: show(rd[i]) if rd.get(i) else 'no decision' for n, rd in zip(names, reads)},
                           'final': show(ad[i]) if i in ad else 'no decision: left as it was',
                           'decision': ad[i][-1] if i in ad else ''})
    return final, agree, review, True


def reader_names(text, swapped, names):
    """The adjudicator's decision with its 'Reading 1' and 'Reading 2' (the order it was shown the readings in, which
    alternates by record) given as the readers' names, so that the decision reads against the review table."""
    first, second = (names[1], names[0]) if swapped else (names[0], names[1])
    text = re.sub(r'\b[Rr]eadings 1 and 2\b', 'both readers', text)
    return re.sub(r'\b[Rr]eading ([12])\b', lambda m: f'reader {first if m.group(1) == "1" else second}', text)
