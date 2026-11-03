import type { Rect, Edge } from './layout';
export const WORKSPACE = Symbol('explorer-workspace');
export type WorkspaceContext = {
  readonly rectangles: Record<string, Rect>;
  readonly dragging: string | null;
  readonly collapsed: string[];
  register(id: string, title: string): () => void;
  begin(id: string, x: number, y: number, title: string): void;
  drag(x: number, y: number): void;
  drop(): void;
  cancel(): void;
  collapse(id: string): void;
  keyboard(id: string, edge: Edge): void;
};
