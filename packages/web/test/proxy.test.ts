import { test } from 'node:test';
import assert from 'node:assert/strict';
import { buildUpstreamUrl, proxyRequest } from '../src/lib/server/proxy-core.ts';

const base = 'http://backend.test:8085';
test('absolute-looking, encoded, and Unicode proxy paths preserve the configured origin', () => {
  for (const path of [
    '/api/https://attacker.test/',
    '/api//attacker.test/',
    '/app-api/\\attacker.test/',
    '/api/%2f%2fattacker.test/',
    '/api/entities/é',
  ]) {
    const target = buildUpstreamUrl(path, '?release=1', base);
    assert.equal(target.origin, base);
    assert.equal(target.search, '?release=1');
  }
  assert.equal(
    buildUpstreamUrl('/api/entities/search', '?q=A', base).href,
    `${base}/entities/search?q=A`,
  );
  assert.throws(() => buildUpstreamUrl('/api/a', '', 'http://user:secret@backend.test'));
});
test('external redirects are rejected before credentials or bodies can leave the backend', async () => {
  const calls: string[] = [];
  const request = new Request('http://web.test/api/jobs', {
    method: 'POST',
    headers: {
      authorization: 'Bearer secret',
      cookie: 'other=private; omnipath_release=latest',
      connection: 'x-admin-secret',
      'x-admin-secret': 'nominated',
    },
    body: 'private body',
  });
  const response = await proxyRequest(request, new URL(request.url), base, async (input, init) => {
    calls.push(String(input));
    const headers = new Headers(init?.headers);
    assert.equal(headers.get('cookie'), 'omnipath_release=latest');
    assert.equal(headers.get('x-admin-secret'), null);
    assert.equal(headers.get('authorization'), 'Bearer secret');
    assert.equal(init?.redirect, 'manual');
    return new Response(null, {
      status: 307,
      headers: { location: 'https://attacker.test/collect' },
    });
  });
  assert.equal(response.status, 502);
  assert.deepEqual(calls, [`${base}/jobs`]);
});
test('same-origin redirects follow HTTP method rules and stay bounded', async () => {
  const methods: string[] = [];
  const request = new Request('http://web.test/api/jobs', { method: 'POST', body: 'data' });
  const response = await proxyRequest(request, new URL(request.url), base, async (_input, init) => {
    methods.push(init!.method!);
    return methods.length === 1
      ? new Response(null, { status: 303, headers: { location: '/result' } })
      : Response.json({ ok: true });
  });
  assert.deepEqual(methods, ['POST', 'GET']);
  assert.equal(response.status, 200);
});
