/** Coalesce synchronous edits and serialize navigations, preserving every parameter. */
export function createUrlUpdater(
  read: () => URL,
  navigate: (url: URL) => Promise<unknown>,
  optimistic: (url: URL) => void,
) {
  let pending: URL | null = null;
  let running = false;
  async function flush() {
    if (running) return;
    running = true;
    try {
      while (pending) {
        const next = pending;
        pending = null;
        await navigate(next);
      }
    } finally {
      running = false;
    }
  }
  return (patch: Record<string, string | null>) => {
    const next = new URL(pending || read());
    for (const [key, value] of Object.entries(patch)) {
      if (value === null) next.searchParams.delete(key);
      else next.searchParams.set(key, value);
    }
    pending = next;
    optimistic(next);
    queueMicrotask(() => {
      void flush().catch((error) => console.error('Selection navigation failed', error));
    });
  };
}
