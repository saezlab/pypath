/** Fixed-origin API proxy. Kept free of framework globals for transport tests. */
const HOP_HEADERS = [
  'connection',
  'keep-alive',
  'proxy-authenticate',
  'proxy-authorization',
  'te',
  'trailer',
  'transfer-encoding',
  'upgrade',
];
const REQUEST_HEADERS = [
  'accept',
  'accept-language',
  'content-type',
  'authorization',
  'x-admin-secret',
  'x-omnipath-release',
  'range',
  'if-range',
  'if-none-match',
  'if-modified-since',
];

export function buildUpstreamUrl(pathname: string, search: string, serviceUrl: string): URL {
  const upstream = new URL(serviceUrl);
  if (!['http:', 'https:'].includes(upstream.protocol) || upstream.username || upstream.password) {
    throw new Error('API_SERVICE_URL must be an HTTP(S) URL without credentials');
  }
  const path = pathname.replace(/^\/(?:app-api|api)(?:\/|$)/, '').replace(/^\/+/, '');
  // Assign the pathname: user input must never be interpreted as a URL origin.
  upstream.pathname = `${upstream.pathname.replace(/\/$/, '')}/${path}`;
  upstream.search = search;
  upstream.hash = '';
  return upstream;
}

function stripHopHeaders(headers: Headers) {
  const nominated = (headers.get('connection') || '').split(',').map((name) => name.trim());
  for (const name of [...HOP_HEADERS, ...nominated]) if (name) headers.delete(name);
}

export async function proxyRequest(
  request: Request,
  url: URL,
  serviceUrl: string,
  fetcher: typeof fetch = fetch,
): Promise<Response> {
  let destination = buildUpstreamUrl(url.pathname, url.search, serviceUrl);
  const origin = destination.origin;
  const incoming = new Headers(request.headers);
  stripHopHeaders(incoming);
  const headers = new Headers();
  for (const name of REQUEST_HEADERS)
    if (incoming.has(name)) headers.set(name, incoming.get(name)!);
  const cookies = (incoming.get('cookie') || '')
    .split(';')
    .map((cookie) => cookie.trim())
    .filter((cookie) => /^omnipath_release=/.test(cookie));
  if (cookies.length) headers.set('cookie', cookies.join('; '));
  headers.set('accept-encoding', 'identity');
  let method = request.method;
  let payload = method === 'GET' || method === 'HEAD' ? undefined : await request.arrayBuffer();
  let upstream: Response;
  for (let redirects = 0; ; redirects++) {
    upstream = await fetcher(destination, {
      method,
      headers,
      body: payload,
      signal: request.signal,
      redirect: 'manual',
    });
    if (![301, 302, 303, 307, 308].includes(upstream.status) || !upstream.headers.has('location'))
      break;
    const next = new URL(upstream.headers.get('location')!, destination);
    await upstream.body?.cancel();
    if (next.origin !== origin || next.username || next.password || redirects >= 4) {
      return Response.json(
        {
          error:
            'Backend redirect is outside the configured API origin or exceeds the redirect limit',
        },
        { status: 502 },
      );
    }
    if (
      (upstream.status === 303 && method !== 'HEAD') ||
      ([301, 302].includes(upstream.status) && method === 'POST')
    ) {
      method = 'GET';
      payload = undefined;
      headers.delete('content-type');
    }
    destination = next;
  }
  const responseHeaders = new Headers(upstream.headers);
  stripHopHeaders(responseHeaders);
  let body = upstream.body;
  if (responseHeaders.has('content-encoding')) {
    responseHeaders.delete('content-encoding');
    responseHeaders.delete('content-length');
  }
  const acceptsGzip = (request.headers.get('accept-encoding') || '').split(',').some((part) => {
    const [coding, ...params] = part.trim().split(';');
    return coding === 'gzip' && !params.some((param) => /^q=0(?:\.0*)?$/.test(param.trim()));
  });
  if (body && responseHeaders.get('content-type')?.includes('application/json')) {
    const vary = responseHeaders.get('vary');
    responseHeaders.set('vary', vary ? `${vary}, Accept-Encoding` : 'Accept-Encoding');
    if (acceptsGzip) {
      body = body.pipeThrough(new CompressionStream('gzip'));
      responseHeaders.set('content-encoding', 'gzip');
      responseHeaders.delete('content-length');
    }
  }
  return new Response(body, {
    status: upstream.status,
    statusText: upstream.statusText,
    headers: responseHeaders,
  });
}
