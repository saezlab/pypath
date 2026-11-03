// The biological topics of the resource catalog, as blocks that select resources.
import type { Component } from 'svelte';
import {
  Activity,
  Apple,
  Dna,
  FlaskConical,
  HeartPulse,
  Layers,
  Microscope,
  Network,
  Pill,
  Shapes,
  Tag,
} from '@lucide/svelte';
import type { ResourceRecord, ResourceTag } from './types';

const icons: Record<string, Component> = {
  signaling: Activity,
  cell_cell_communication: Network,
  gene_regulation: Dna,
  metabolism: FlaskConical,
  subcellular_localization: Layers,
  disease_phenotype: HeartPulse,
  pharmacology: Pill,
  molecular_classification: Shapes,
  cell_types: Microscope,
  nutrition: Apple,
};

export function topicIcon(id: string): Component {
  return icons[id] ?? Tag;
}

export interface Topic {
  tag: ResourceTag;
  /** The resources with the topic, most relations first. */
  members: ResourceRecord[];
}

/** One topic per `topic` tag of the catalog, the topics with most resources first. */
export function catalogTopics(resources: ResourceRecord[]): Topic[] {
  const byTopic = new Map<string, { tag: ResourceTag; members: ResourceRecord[] }>();
  for (const resource of resources)
    for (const tag of resource.tags ?? [])
      if (tag.dimension === 'topic') {
        const entry = byTopic.get(tag.id) ?? { tag, members: [] };
        entry.members.push(resource);
        byTopic.set(tag.id, entry);
      }
  return [...byTopic.values()]
    .map(({ tag, members }) => ({
      tag,
      members: [...members].sort((a, b) => b.interaction_count - a.interaction_count),
    }))
    .sort((a, b) => b.members.length - a.members.length || a.tag.label.localeCompare(b.tag.label));
}
