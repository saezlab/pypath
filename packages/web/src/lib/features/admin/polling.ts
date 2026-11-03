/** Schedule the next poll after completion; manual refreshes share ongoing work. */
export function createPoller(
  task: (signal: AbortSignal) => Promise<unknown>,
  delay: () => number,
  onError: (error: unknown) => void = console.error,
) {
  let active = false;
  let generation = 0;
  let timer: ReturnType<typeof setTimeout> | undefined;
  let controller: AbortController | undefined;
  let pending: Promise<unknown> | undefined;
  function refresh(): Promise<unknown> {
    if (pending) return pending;
    controller = new AbortController();
    const signal = controller.signal;
    pending = Promise.resolve()
      .then(() => task(signal))
      .catch((error) => {
        if (!signal.aborted) onError(error);
      })
      .finally(() => {
        pending = undefined;
      });
    return pending;
  }
  async function tick(current: number) {
    await refresh();
    if (active && current === generation)
      timer = setTimeout(() => {
        void tick(current);
      }, delay());
  }
  return {
    refresh,
    start() {
      if (!active) {
        active = true;
        void tick(++generation);
      }
    },
    stop() {
      active = false;
      generation++;
      clearTimeout(timer);
      controller?.abort();
    },
  };
}
