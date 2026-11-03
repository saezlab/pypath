declare module 'canvas-confetti' {
  interface Options {
    particleCount?: number;
    angle?: number;
    spread?: number;
    startVelocity?: number;
    decay?: number;
    gravity?: number;
    drift?: number;
    ticks?: number;
    origin?: {
      x?: number;
      y?: number;
    };
    colors?: string[];
    shapes?: ('square' | 'circle' | 'star')[];
    zIndex?: number;
    disableForReducedMotion?: boolean;
    scalar?: number;
  }

  interface ConfettiFunction {
    (options?: Options): Promise<null>;
    reset(): void;
  }

  const confetti: ConfettiFunction;
  export default confetti;
}
