/** Shows the full text as a native tooltip, but only while the element's text is cut off. */
export function truncateTitle(node: HTMLElement, text?: string | null) {
  let full = text ?? '';
  const update = () => {
    const truncated =
      node.scrollWidth > node.clientWidth + 1 || node.scrollHeight > node.clientHeight + 1;
    const value = full || node.textContent?.trim() || '';
    if (truncated && value) node.title = value;
    else node.removeAttribute('title');
  };
  node.addEventListener('pointerenter', update);
  return {
    update(next?: string | null) {
      full = next ?? '';
    },
    destroy() {
      node.removeEventListener('pointerenter', update);
    },
  };
}
