<script lang="ts">
  import { getIdentifierTypeLabel } from '$lib/domain/display';

  interface Props {
    identifierType: string;
    value: string;
    variant?: 'default' | 'subtle' | 'compact';
    class?: string;
  }

  let { identifierType, value, variant = 'default', class: className = '' }: Props = $props();

  const typeLabel = $derived(getIdentifierTypeLabel(identifierType));
  const rootClass = $derived(
    variant === 'compact'
      ? 'inline-flex h-5 max-w-full items-center overflow-hidden rounded-md border border-border/70 bg-background/80 text-[10px] leading-none shadow-xs'
      : variant === 'subtle'
        ? 'inline-flex h-6 max-w-full items-center overflow-hidden rounded-md border border-border/70 bg-background/80 text-xs shadow-xs'
        : 'inline-flex h-7 max-w-full items-center overflow-hidden rounded-md border border-border/80 bg-card text-xs shadow-xs',
  );
  const typeClass = $derived(
    variant === 'compact'
      ? 'flex h-full max-w-16 shrink-0 items-center truncate border-r border-border/60 bg-muted/45 px-1.5 font-medium text-muted-foreground'
      : variant === 'subtle'
        ? 'flex h-full max-w-24 shrink-0 items-center truncate border-r border-border/60 bg-muted/45 px-2 font-medium text-muted-foreground'
        : 'flex h-full max-w-28 shrink-0 items-center truncate border-r border-border/70 bg-muted/50 px-2.5 font-medium text-muted-foreground',
  );
  const valueClass = $derived(
    variant === 'compact'
      ? 'min-w-0 truncate px-1.5 font-mono text-foreground'
      : variant === 'subtle'
        ? 'min-w-0 truncate px-2 font-mono text-foreground'
        : 'min-w-0 truncate px-2.5 font-mono text-foreground',
  );
</script>

<span class={`${rootClass} ${className}`} title={`${typeLabel}: ${value}`}>
  <span class={typeClass}>
    {typeLabel}
  </span>
  <span class={valueClass}>
    {value}
  </span>
</span>
