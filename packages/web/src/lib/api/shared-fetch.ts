/** Coalesce in-flight reads; each consumer can cancel without cancelling others. */
const pending = new Map<
  string,
  { promise: Promise<Response>; controller: AbortController; users: number }
>();

export function sharedFetch(input: string | URL, init: RequestInit = {}): Promise<Response> {
  const key = JSON.stringify([
    String(input),
    init.method || 'GET',
    [...new Headers(init.headers)],
    init.body || null,
  ]);
  const signal = init.signal;
  if (signal?.aborted) return Promise.reject(new DOMException('Aborted', 'AbortError'));
  let entry = pending.get(key);
  if (!entry) {
    const controller = new AbortController();
    entry = { controller, users: 0, promise: fetch(input, { ...init, signal: controller.signal }) };
    pending.set(key, entry);
    const current = entry;
    const cleanup = () => {
      if (pending.get(key) === current) pending.delete(key);
    };
    entry.promise.then(cleanup, cleanup);
  }
  entry.users++;
  const current = entry;
  return new Promise((resolve, reject) => {
    let finished = false;
    const release = () => {
      if (finished) return false;
      finished = true;
      signal?.removeEventListener('abort', abort);
      current.users--;
      if (!current.users)
        queueMicrotask(() => {
          if (!current.users && pending.get(key) === current) {
            pending.delete(key);
            current.controller.abort();
          }
        });
      return true;
    };
    const abort = () => {
      if (release()) reject(new DOMException('Aborted', 'AbortError'));
    };
    signal?.addEventListener('abort', abort, { once: true });
    current.promise.then(
      (response) => {
        if (release()) resolve(response.clone());
      },
      (error) => {
        if (release()) reject(error);
      },
    );
  });
}
