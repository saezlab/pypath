<script lang="ts">
  import { onDestroy } from 'svelte';
  import { page } from '$app/state';
  import { Check, Copy, Database, Download, GitPullRequest } from '@lucide/svelte';
  import { toast } from 'svelte-sonner';
  import { Button } from '$lib/components/ui/button/index.js';
  import { skills } from '$lib/skills';

  // Skills are listed in the order most visitors need them; the first is where to start.
  let selectedName = $state(skills[0].name);
  const selected = $derived(skills.find((skill) => skill.name === selectedName) ?? skills[0]);

  // A prompt any agent can act on: it fetches the hosted SKILL.md and installs it where it
  // loads skills from.
  const prompt = $derived(
    `Set up the "${selected.title}" agent skill from ${page.url.origin}${selected.href}. ` +
      `Save it as ${selected.name}/SKILL.md where you load skills from (for Claude Code: ~/.claude/skills), ` +
      'or, if you do not use skills, in this workspace and refer to it from AGENTS.md.',
  );

  let copied = $state(false);
  let copyTimer: ReturnType<typeof setTimeout> | undefined;
  onDestroy(() => clearTimeout(copyTimer));

  async function copyPrompt() {
    try {
      await navigator.clipboard.writeText(prompt);
      copied = true;
      toast.success('Prompt copied. Paste it into your agent.');
      clearTimeout(copyTimer);
      copyTimer = setTimeout(() => (copied = false), 1800);
    } catch {
      toast.error('Could not copy. Select the prompt and copy it instead.');
    }
  }
</script>

<svelte:head>
  <title>Skills · OmniPath</title>
  <meta
    name="description"
    content="Teach your agent to use OmniPath data or contribute a new resource to pypath. Copy a setup prompt or download a skill."
  />
</svelte:head>

{#snippet step(number: number, title: string)}
  <span
    class="flex size-6 items-center justify-center rounded-full bg-foreground text-xs font-semibold text-background"
    aria-hidden="true">{number}</span
  >
  <h2 class="text-[15px] font-semibold">{title}</h2>
{/snippet}

<div class="max-w-4xl space-y-8 pt-4 pb-10">
  <h1 class="text-2xl font-semibold tracking-tight">Agent skills</h1>

  <section class="grid grid-cols-[1.5rem_minmax(0,1fr)] items-center gap-3">
    {@render step(1, 'Choose a skill')}
    <fieldset class="col-start-2 grid gap-2.5 sm:grid-cols-2">
      <legend class="sr-only">Skill</legend>
      {#each skills as skill, index (skill.name)}
        <label
          class="grid cursor-pointer grid-cols-[auto_minmax(0,1fr)] gap-x-3 gap-y-1 rounded-xl border bg-background p-3.5 transition-colors hover:border-primary/50 has-checked:border-primary has-checked:bg-primary/5 has-checked:ring-1 has-checked:ring-primary has-focus-visible:ring-2 has-focus-visible:ring-ring"
        >
          <input
            type="radio"
            name="skill"
            value={skill.name}
            bind:group={selectedName}
            class="sr-only"
          />
          <span
            class="row-span-2 flex size-9 items-center justify-center rounded-lg bg-primary/10 text-primary"
            aria-hidden="true"
          >
            {#if skill.kind === 'Contribute'}<GitPullRequest class="size-4" />{:else}<Database
                class="size-4"
              />{/if}
          </span>
          <span class="flex flex-wrap items-center gap-2 text-[15px] font-semibold">
            {skill.title}
            {#if index === 0}
              <span
                class="rounded-full bg-primary/10 px-2 py-0.5 text-[11px] font-semibold text-primary"
                >Start here</span
              >
            {/if}
          </span>
          <span class="text-sm text-muted-foreground">{skill.description}</span>
        </label>
      {/each}
    </fieldset>
  </section>

  <section class="grid grid-cols-[1.5rem_minmax(0,1fr)] items-center gap-3">
    {@render step(2, 'Set it up')}
    <div class="col-start-2 space-y-3">
      <p class="rounded-lg border bg-muted/30 px-3.5 py-3 text-xs leading-relaxed">{prompt}</p>
      <div class="flex flex-wrap gap-2">
        <Button size="sm" onclick={copyPrompt}>
          {#if copied}<Check />Copied{:else}<Copy />Copy prompt for your agent{/if}
        </Button>
        <Button
          size="sm"
          variant="outline"
          href={selected.href}
          download={`${selected.name}-SKILL.md`}><Download />Download SKILL.md</Button
        >
      </div>
    </div>
  </section>
</div>
