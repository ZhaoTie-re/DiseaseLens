#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The language-model calls of the pipeline, isolated and cached.
#           Each call is `claude -p` with every customisation off (--safe-mode:
#           no CLAUDE.md, skills, plugins, hooks, MCP servers), no tools, no
#           session saved, auto-memory disabled, run from a scratch directory,
#           so no project memory or file reaches the model: it sees only the
#           prompt. The answer must be one JSON value.
#           Cache: one file per (model, effort, system prompt, prompt) hash under
#           --cache-dir; a cached answer is reused, so a later snapshot calls the
#           model only for new or changed inputs. Offline (cache only): the model
#           is never called; a question without a cached answer gets None.
#           Grounding helpers check an answer against the text it came from.
# ---------------------------------------------------------------------------
import hashlib
import json
import os
import re
import subprocess
import time
import unicodedata
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


class LLM:
    def __init__(self, cache_dir, workdir, model='', effort='', timeout=600, parallel=3, claude='claude',
                 offline=False):
        self.cache = Path(cache_dir)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.model, self.effort, self.timeout, self.parallel, self.claude = model, effort, timeout, parallel, claude
        self.offline = offline
        self.calls, self.cached, self.missed = 0, 0, 0

    def key(self, system, prompt):
        h = hashlib.sha256('\x1e'.join([self.model, self.effort, system, prompt]).encode()).hexdigest()
        return self.cache / h[:2] / f'{h}.json'

    def ask(self, system, prompt, tries=3):
        """One answer (a parsed JSON value), from the cache or the model."""
        f = self.key(system, prompt)
        if f.exists():
            self.cached += 1
            return json.loads(f.read_text())['answer']
        if self.offline:
            self.missed += 1
            return None
        cmd = [self.claude, '-p', '--safe-mode', '--tools', '', '--strict-mcp-config', '--no-session-persistence',
               '--output-format', 'json', '--system-prompt', system]
        if self.model:
            cmd += ['--model', self.model]
        if self.effort:
            cmd += ['--effort', self.effort]
        env = {**os.environ, 'CLAUDE_CODE_DISABLE_AUTO_MEMORY': '1'}
        env.pop('CLAUDECODE', None)
        err = ''
        for k in range(tries):
            r = subprocess.run(cmd, input=prompt, capture_output=True, text=True, cwd=self.workdir, env=env,
                               timeout=self.timeout)
            try:
                out = json.loads(r.stdout)
                if out.get('is_error'):
                    raise ValueError(out.get('result', '')[:300])
                answer = parse_json(out['result'])
                f.parent.mkdir(parents=True, exist_ok=True)
                f.write_text(json.dumps({'answer': answer, 'model': list((out.get('modelUsage') or {}).keys()),
                                         'cost_usd': out.get('total_cost_usd'), 'at': time.strftime('%F %T')},
                                        ensure_ascii=False))
                self.calls += 1
                return answer
            except (ValueError, KeyError, json.JSONDecodeError) as e:
                err = f'{e}; stderr: {r.stderr[-300:]}'
                time.sleep(5 * (k + 1))
        raise RuntimeError(f'no JSON answer after {tries} tries: {err}')

    def map(self, jobs):
        """[(system, prompt), ...] -> answers in order, `parallel` calls at a time."""
        with ThreadPoolExecutor(self.parallel) as ex:
            return list(ex.map(lambda j: self.ask(*j), jobs))


def parse_json(text):
    """The JSON value of an answer, with or without a ```json fence."""
    t = text.strip()
    m = re.search(r'```(?:json)?\s*(.*?)```', t, re.S)
    if m:
        t = m.group(1).strip()
    if not t.startswith(('{', '[')):
        i = min([x for x in (t.find('{'), t.find('[')) if x >= 0], default=-1)
        if i < 0:
            raise ValueError('no JSON in the answer')
        t = t[i:]
    return json.loads(t)


# ── grounding ────────────────────────────────────────────────────────────────
def norm(s):
    """Text for substring checks: Unicode-normalised, dashes and spaces unified, lower case."""
    s = unicodedata.normalize('NFKC', str(s or ''))
    s = re.sub(r'[‐-―−]', '-', s)
    s = s.replace(' ', ' ').replace('’', "'").replace('“', '"').replace('”', '"')
    return ' '.join(s.split()).lower()


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
