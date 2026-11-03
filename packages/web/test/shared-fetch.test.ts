import { test } from 'node:test';
import assert from 'node:assert/strict';
import { sharedFetch } from '../src/lib/api/shared-fetch';

test('identical reads share transport while cancellation is consumer-specific', async () => {
  const original = globalThis.fetch;
  let finish!: (value: Response) => void;
  let calls = 0;
  let transportSignal: AbortSignal | undefined | null;
  globalThis.fetch = ((_input, init) => {
    calls++;
    transportSignal = init?.signal;
    return new Promise<Response>((resolve) => {
      finish = resolve;
    });
  }) as typeof fetch;
  try {
    const a = new AbortController();
    const first = sharedFetch('/same', { signal: a.signal });
    const second = sharedFetch('/same');
    a.abort();
    await assert.rejects(first, { name: 'AbortError' });
    assert.equal(transportSignal?.aborted, false);
    finish(new Response('{"ok":true}'));
    assert.deepEqual(await (await second).json(), { ok: true });
    assert.equal(calls, 1);
  } finally {
    globalThis.fetch = original;
  }
});

test('release header separates requests and abandoned transport is aborted', async () => {
  const original = globalThis.fetch;
  const signals: AbortSignal[] = [];
  globalThis.fetch = ((_input, init) => {
    signals.push(init!.signal!);
    return new Promise<Response>((_, reject) =>
      init!.signal!.addEventListener('abort', () =>
        reject(new DOMException('Aborted', 'AbortError')),
      ),
    );
  }) as typeof fetch;
  try {
    const a = new AbortController(),
      b = new AbortController();
    const first = sharedFetch('/same', {
      headers: { 'X-OmniPath-Release': '1' },
      signal: a.signal,
    });
    const second = sharedFetch('/same', {
      headers: { 'X-OmniPath-Release': '2' },
      signal: b.signal,
    });
    a.abort();
    b.abort();
    await Promise.all([assert.rejects(first), assert.rejects(second)]);
    assert.equal(signals.length, 2);
    assert.ok(signals.every((signal) => signal.aborted));
  } finally {
    globalThis.fetch = original;
  }
});
