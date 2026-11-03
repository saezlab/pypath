<script lang="ts">
  import { FlaskConical } from '@lucide/svelte';
  import { onMount, tick } from 'svelte';
  import { Alert, AlertDescription } from '$lib/components/ui/alert/index.js';
  import { Button } from '$lib/components/ui/button/index.js';
  import { cn } from '$lib/utils';

  type OpenChemLibModule = {
    Molecule: {
      fromSmiles(smiles: string): {
        toSVG(width: number, height: number): string;
      } | null;
    };
  };

  interface Props {
    smiles: string;
    width?: number;
    height?: number;
    class?: string;
    renderOnClick?: boolean;
  }

  let {
    smiles,
    width = 180,
    height = 140,
    class: className = '',
    renderOnClick = true,
  }: Props = $props();

  let container: HTMLDivElement | null = $state(null);
  let error = $state<string | null>(null);
  let shouldRender = $state(false);
  let hasRendered = $state(false);
  let renderedSmiles = $state<string | null>(null);
  let isVisible = $state(false);
  let ocl = $state<OpenChemLibModule | null>(null);
  let observer: IntersectionObserver | null = null;

  const dimensions = $derived(`width: ${width}px; height: ${height}px;`);

  function isBlackColor(value: string | null) {
    if (!value) return false;
    const normalized = value.trim().toLowerCase();
    return [
      'black',
      '#000',
      '#000000',
      'rgb(0,0,0)',
      'rgb(0, 0, 0)',
      'rgba(0,0,0,1)',
      'rgba(0, 0, 0, 1)',
    ].includes(normalized);
  }

  async function loadOpenChemLib() {
    if (ocl) return;
    try {
      error = null;
      const { Molecule } = await import('openchemlib');
      ocl = { Molecule };
    } catch {
      error = 'Failed to load molecular visualization library';
    }
  }

  function applyDarkModeColors(svgElement: SVGSVGElement) {
    if (!document.documentElement.classList.contains('dark')) return;

    for (const node of svgElement.querySelectorAll('*')) {
      const stroke = node.getAttribute('stroke');
      const fill = node.getAttribute('fill');

      if (isBlackColor(stroke)) node.setAttribute('stroke', '#E5E7EB');
      if (isBlackColor(fill)) node.setAttribute('fill', '#E5E7EB');
    }
  }

  function renderMolecule() {
    if (!container || !ocl || !smiles || hasRendered) return;

    try {
      error = null;
      container.innerHTML = '';

      const molecule = ocl.Molecule.fromSmiles(smiles);
      if (!molecule) {
        error = 'Invalid molecular structure';
        return;
      }

      const svgString = molecule.toSVG(width, height);
      if (!svgString) {
        error = 'Failed to generate structure visualization';
        return;
      }

      container.innerHTML = svgString;
      const svgElement = container.querySelector('svg');
      if (svgElement) {
        svgElement.style.width = '100%';
        svgElement.style.height = '100%';
        svgElement.style.display = 'block';
        applyDarkModeColors(svgElement);
      }
      hasRendered = true;
      renderedSmiles = smiles;
    } catch {
      error = 'Failed to render molecular structure';
    }
  }

  onMount(() => {
    return () => {
      observer?.disconnect();
    };
  });

  $effect(() => {
    if (!renderOnClick) shouldRender = true;
  });

  $effect(() => {
    if (!renderedSmiles || renderedSmiles === smiles) return;

    hasRendered = false;
    renderedSmiles = null;
    if (container) container.innerHTML = '';
  });

  $effect(() => {
    if (!container || !shouldRender || hasRendered) return;

    observer?.disconnect();
    observer = new IntersectionObserver(
      (entries) => {
        if (!entries[0]?.isIntersecting) return;
        isVisible = true;
        observer?.disconnect();
      },
      { rootMargin: '300px 0px' },
    );
    observer.observe(container);

    return () => observer?.disconnect();
  });

  $effect(() => {
    if (!shouldRender || !isVisible) return;

    loadOpenChemLib();
  });

  $effect(() => {
    if (!shouldRender || !isVisible || !ocl || hasRendered) return;

    tick().then(renderMolecule);
  });
</script>

<div class="flex flex-col items-center" style={`width: ${width}px;`}>
  {#if error}
    <Alert class={cn('flex items-center justify-center text-center', className)} style={dimensions}>
      <AlertDescription class="text-sm">{error}</AlertDescription>
    </Alert>
  {:else}
    <div
      bind:this={container}
      class={cn(
        'relative shrink-0 overflow-hidden rounded-md border bg-muted/20',
        !hasRendered && shouldRender ? 'animate-pulse' : '',
        className,
      )}
      style={dimensions}
    >
      {#if !shouldRender}
        <div class="absolute inset-0 flex items-center justify-center bg-muted/10">
          <Button
            type="button"
            variant="secondary"
            size="sm"
            class="rounded-full"
            onclick={() => {
              shouldRender = true;
            }}
          >
            <FlaskConical class="mr-1.5 size-4" />
            Show structure
          </Button>
        </div>
      {/if}
    </div>
  {/if}
</div>
