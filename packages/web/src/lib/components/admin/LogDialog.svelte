<script lang="ts">
  import { Check, Copy } from '@lucide/svelte';
  import { Button } from '$lib/components/ui/button/index.js';
  import {
    Dialog,
    DialogContent,
    DialogHeader,
    DialogTitle,
  } from '$lib/components/ui/dialog/index.js';

  interface Props {
    open?: boolean;
    title: string;
    text: string;
    loading?: boolean;
  }

  let { open = $bindable(false), title, text, loading = false }: Props = $props();
  let copied = $state(false);
  let copyTimer: ReturnType<typeof setTimeout> | null = null;

  async function copyLogs() {
    try {
      await navigator.clipboard.writeText(text || '');
      copied = true;
      if (copyTimer) clearTimeout(copyTimer);
      copyTimer = setTimeout(() => {
        copied = false;
      }, 1600);
    } catch {
      copied = false;
    }
  }
</script>

<Dialog bind:open>
  <DialogContent
    class="grid max-h-[85vh] grid-rows-[auto_minmax(0,1fr)] gap-3 overflow-hidden sm:max-w-3xl"
  >
    <DialogHeader class="pr-10">
      <div class="flex items-center justify-between gap-3">
        <DialogTitle class="truncate">{title}</DialogTitle>
        <Button size="sm" variant="outline" onclick={copyLogs} disabled={loading || !text}>
          {#if copied}
            <Check class="size-3.5" />
            Copied
          {:else}
            <Copy class="size-3.5" />
            Copy
          {/if}
        </Button>
      </div>
    </DialogHeader>
    <div class="min-h-0 overflow-auto rounded-lg bg-zinc-950 p-4 text-zinc-100">
      {#if loading}
        <p class="text-xs text-zinc-400">Loading logs…</p>
      {:else if text}
        <pre class="whitespace-pre-wrap break-words font-mono text-xs leading-relaxed">{text}</pre>
      {:else}
        <p class="text-xs text-zinc-400">No logs saved yet.</p>
      {/if}
    </div>
  </DialogContent>
</Dialog>
