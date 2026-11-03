let width = $state(1280);
export const workspaceWidth = {
  get value() {
    return width;
  },
  load() {
    try {
      const saved = Number(localStorage.getItem('omnipath-workspace-width'));
      if (saved >= 960 && saved <= 3840) width = saved;
    } catch {}
  },
  set(value: number) {
    width = Math.max(960, Math.min(3840, value));
    try {
      localStorage.setItem('omnipath-workspace-width', String(width));
    } catch {}
  },
};
