export interface PageState<T> {
  items: T[];
  loading: boolean;
  loadingMore: boolean;
  hasMore: boolean;
  error: string | null;
}
export type PageLoader<T> = (offset: number, signal: AbortSignal) => Promise<T[]>;

/** A query owns every page and abort signal; replaced queries cannot publish. */
export function createPagedQuery<T>(pageSize: number, publish: (state: PageState<T>) => void) {
  let generation = 0;
  let controller: AbortController | null = null;
  let loader: PageLoader<T> | null = null;
  let state: PageState<T> = {
    items: [],
    loading: false,
    loadingMore: false,
    hasMore: true,
    error: null,
  };
  const update = (patch: Partial<PageState<T>>) => {
    state = { ...state, ...patch };
    publish(state);
  };
  async function load(first: boolean) {
    if (
      !loader ||
      !controller ||
      (!first && (state.loading || state.loadingMore || !state.hasMore))
    )
      return;
    const current = generation;
    const signal = controller.signal;
    const offset = first ? 0 : state.items.length;
    update(first ? { loading: true, error: null } : { loadingMore: true, error: null });
    try {
      const page = await loader(offset, signal);
      if (signal.aborted || current !== generation) return;
      update({
        items: first ? page : [...state.items, ...page],
        hasMore: page.length === pageSize,
      });
    } catch (error) {
      if (signal.aborted || current !== generation) return;
      update({
        error: error instanceof Error ? error.message : 'Failed to load results',
        ...(first ? { hasMore: false } : {}),
      });
    } finally {
      if (current === generation && !signal.aborted) update({ loading: false, loadingMore: false });
    }
  }
  function cancel() {
    generation++;
    controller?.abort();
    controller = null;
  }
  return {
    reset(next: PageLoader<T>) {
      cancel();
      controller = new AbortController();
      loader = next;
      update({ items: [], loading: true, loadingMore: false, hasMore: true, error: null });
      return load(true);
    },
    loadMore: () => load(false),
    cancel,
  };
}
