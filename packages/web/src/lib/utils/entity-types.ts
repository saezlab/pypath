import { BIOLINK_ENTITIES, ENTITY_LOOKUP } from '$lib/domain/biolink.generated';

// Presentation choices keyed by canonical Biolink class names.
export const entityTypeEmojis: Record<string, string> = {
  anatomical_entity: '🫀',
  cell: '🦠',
  food: '🍎',
  gene_family: '👥',
  information_content_entity: '📄',
  microrna: '🧵',
  molecular_entity: '⚛️',
  organism_taxon: '🌿',
  physical_entity: '📦',
  polypeptide: '🔗',
  chemical_entity: '🧪',
  small_molecule: '🧪',
  ontology_class: '🏷️',
  protein: '🧬',
  protein_family: '👥',
  gene: '🧬',
  rna_product: '🧫',
  nucleic_acid_entity: '🧬',
  macromolecular_complex: '🧩',
  pathway: '🛣️',
  biological_process: '⚙️',
  molecular_activity: '⚗️',
  cellular_component: '🧱',
  phenotypic_feature: '🩺',
  exposure_event: '🔦',
  named_thing: '🔹',
};

export type EntityTypeStyle = {
  label: string;
  color: string;
  bgColor: string;
  borderColor: string;
  chipClass: string;
  /** Faint tint for a result row, with a stronger hover. */
  rowClass: string;
  /** Icon tile on a result row. */
  iconClass: string;
};

export const defaultEntityTypeStyle: EntityTypeStyle = {
  label: 'Entity',
  color: 'text-slate-500',
  bgColor: 'from-slate-50/80 to-slate-100/80 dark:from-slate-800/80 dark:to-slate-900/80',
  borderColor: 'border-slate-200/60 dark:border-slate-700/60',
  chipClass:
    'border-slate-200 bg-slate-50 text-slate-700 dark:border-slate-700 dark:bg-slate-900/45 dark:text-slate-200',
  rowClass:
    'bg-slate-500/[0.05] hover:bg-slate-500/[0.10] dark:bg-slate-400/[0.06] dark:hover:bg-slate-400/[0.11]',
  iconClass: 'bg-slate-500/15 text-slate-600 dark:bg-slate-400/15 dark:text-slate-300',
};

const entityTypeStyles: Record<string, EntityTypeStyle> = {
  protein: {
    label: 'Protein',
    color: 'text-blue-500',
    bgColor: 'from-blue-50/80 to-blue-100/80 dark:from-blue-900/30 dark:to-blue-800/30',
    borderColor: 'border-blue-200/70 dark:border-blue-700/60',
    chipClass:
      'border-blue-200 bg-blue-50 text-blue-700 dark:border-blue-700/60 dark:bg-blue-950/35 dark:text-blue-200',
    rowClass:
      'bg-blue-500/[0.05] hover:bg-blue-500/[0.10] dark:bg-blue-400/[0.06] dark:hover:bg-blue-400/[0.11]',
    iconClass: 'bg-blue-500/15 text-blue-600 dark:bg-blue-400/15 dark:text-blue-300',
  },
  chemical_entity: {
    label: 'Chemical entity',
    color: 'text-green-500',
    bgColor: 'from-green-50/80 to-green-100/80 dark:from-green-900/30 dark:to-green-800/30',
    borderColor: 'border-green-200/70 dark:border-green-700/60',
    chipClass:
      'border-green-200 bg-green-50 text-green-700 dark:border-green-700/60 dark:bg-green-950/35 dark:text-green-200',
    rowClass:
      'bg-green-500/[0.05] hover:bg-green-500/[0.10] dark:bg-green-400/[0.06] dark:hover:bg-green-400/[0.11]',
    iconClass: 'bg-green-500/15 text-green-600 dark:bg-green-400/15 dark:text-green-300',
  },
  compound: {
    label: 'Compound',
    color: 'text-green-500',
    bgColor: 'from-green-50/80 to-green-100/80 dark:from-green-900/30 dark:to-green-800/30',
    borderColor: 'border-green-200/70 dark:border-green-700/60',
    chipClass:
      'border-green-200 bg-green-50 text-green-700 dark:border-green-700/60 dark:bg-green-950/35 dark:text-green-200',
    rowClass:
      'bg-green-500/[0.05] hover:bg-green-500/[0.10] dark:bg-green-400/[0.06] dark:hover:bg-green-400/[0.11]',
    iconClass: 'bg-green-500/15 text-green-600 dark:bg-green-400/15 dark:text-green-300',
  },
  metabolite: {
    label: 'Metabolite',
    color: 'text-green-500',
    bgColor: 'from-green-50/80 to-green-100/80 dark:from-green-900/30 dark:to-green-800/30',
    borderColor: 'border-green-200/70 dark:border-green-700/60',
    chipClass:
      'border-green-200 bg-green-50 text-green-700 dark:border-green-700/60 dark:bg-green-950/35 dark:text-green-200',
    rowClass:
      'bg-green-500/[0.05] hover:bg-green-500/[0.10] dark:bg-green-400/[0.06] dark:hover:bg-green-400/[0.11]',
    iconClass: 'bg-green-500/15 text-green-600 dark:bg-green-400/15 dark:text-green-300',
  },
  drug: {
    label: 'Drug',
    color: 'text-purple-500',
    bgColor: 'from-purple-50/80 to-purple-100/80 dark:from-purple-900/30 dark:to-purple-800/30',
    borderColor: 'border-purple-200/70 dark:border-purple-700/60',
    chipClass:
      'border-purple-200 bg-purple-50 text-purple-700 dark:border-purple-700/60 dark:bg-purple-950/35 dark:text-purple-200',
    rowClass:
      'bg-purple-500/[0.05] hover:bg-purple-500/[0.10] dark:bg-purple-400/[0.06] dark:hover:bg-purple-400/[0.11]',
    iconClass: 'bg-purple-500/15 text-purple-600 dark:bg-purple-400/15 dark:text-purple-300',
  },
  lipid: {
    label: 'Lipid',
    color: 'text-yellow-600',
    bgColor: 'from-yellow-50/80 to-yellow-100/80 dark:from-yellow-900/30 dark:to-yellow-800/30',
    borderColor: 'border-yellow-200/80 dark:border-yellow-700/60',
    chipClass:
      'border-yellow-200 bg-yellow-50 text-yellow-800 dark:border-yellow-700/60 dark:bg-yellow-950/35 dark:text-yellow-200',
    rowClass:
      'bg-yellow-500/[0.05] hover:bg-yellow-500/[0.10] dark:bg-yellow-400/[0.06] dark:hover:bg-yellow-400/[0.11]',
    iconClass: 'bg-yellow-500/15 text-yellow-600 dark:bg-yellow-400/15 dark:text-yellow-300',
  },
  gene: {
    label: 'Gene',
    color: 'text-orange-500',
    bgColor: 'from-orange-50/80 to-orange-100/80 dark:from-orange-900/30 dark:to-orange-800/30',
    borderColor: 'border-orange-200/70 dark:border-orange-700/60',
    chipClass:
      'border-orange-200 bg-orange-50 text-orange-700 dark:border-orange-700/60 dark:bg-orange-950/35 dark:text-orange-200',
    rowClass:
      'bg-orange-500/[0.05] hover:bg-orange-500/[0.10] dark:bg-orange-400/[0.06] dark:hover:bg-orange-400/[0.11]',
    iconClass: 'bg-orange-500/15 text-orange-600 dark:bg-orange-400/15 dark:text-orange-300',
  },
  macromolecular_complex: {
    label: 'Macromolecular complex',
    color: 'text-indigo-500',
    bgColor: 'from-indigo-50/80 to-indigo-100/80 dark:from-indigo-900/30 dark:to-indigo-800/30',
    borderColor: 'border-indigo-200/70 dark:border-indigo-700/60',
    chipClass:
      'border-indigo-200 bg-indigo-50 text-indigo-700 dark:border-indigo-700/60 dark:bg-indigo-950/35 dark:text-indigo-200',
    rowClass:
      'bg-indigo-500/[0.05] hover:bg-indigo-500/[0.10] dark:bg-indigo-400/[0.06] dark:hover:bg-indigo-400/[0.11]',
    iconClass: 'bg-indigo-500/15 text-indigo-600 dark:bg-indigo-400/15 dark:text-indigo-300',
  },
  pathway: {
    label: 'Pathway',
    color: 'text-cyan-500',
    bgColor: 'from-cyan-50/80 to-cyan-100/80 dark:from-cyan-900/30 dark:to-cyan-800/30',
    borderColor: 'border-cyan-200/70 dark:border-cyan-700/60',
    chipClass:
      'border-cyan-200 bg-cyan-50 text-cyan-700 dark:border-cyan-700/60 dark:bg-cyan-950/35 dark:text-cyan-200',
    rowClass:
      'bg-cyan-500/[0.05] hover:bg-cyan-500/[0.10] dark:bg-cyan-400/[0.06] dark:hover:bg-cyan-400/[0.11]',
    iconClass: 'bg-cyan-500/15 text-cyan-600 dark:bg-cyan-400/15 dark:text-cyan-300',
  },
  reaction: {
    label: 'Reaction',
    color: 'text-pink-500',
    bgColor: 'from-pink-50/80 to-pink-100/80 dark:from-pink-900/30 dark:to-pink-800/30',
    borderColor: 'border-pink-200/70 dark:border-pink-700/60',
    chipClass:
      'border-pink-200 bg-pink-50 text-pink-700 dark:border-pink-700/60 dark:bg-pink-950/35 dark:text-pink-200',
    rowClass:
      'bg-pink-500/[0.05] hover:bg-pink-500/[0.10] dark:bg-pink-400/[0.06] dark:hover:bg-pink-400/[0.11]',
    iconClass: 'bg-pink-500/15 text-pink-600 dark:bg-pink-400/15 dark:text-pink-300',
  },
  molecular_activity: {
    label: 'Molecular activity',
    color: 'text-pink-500',
    bgColor: 'from-pink-50/80 to-pink-100/80 dark:from-pink-900/30 dark:to-pink-800/30',
    borderColor: 'border-pink-200/70 dark:border-pink-700/60',
    chipClass:
      'border-pink-200 bg-pink-50 text-pink-700 dark:border-pink-700/60 dark:bg-pink-950/35 dark:text-pink-200',
    rowClass:
      'bg-pink-500/[0.05] hover:bg-pink-500/[0.10] dark:bg-pink-400/[0.06] dark:hover:bg-pink-400/[0.11]',
    iconClass: 'bg-pink-500/15 text-pink-600 dark:bg-pink-400/15 dark:text-pink-300',
  },
  ontology_class: {
    label: 'Ontology class',
    color: 'text-amber-600',
    bgColor: 'from-amber-50/80 to-amber-100/80 dark:from-amber-900/30 dark:to-amber-800/30',
    borderColor: 'border-amber-200/80 dark:border-amber-700/60',
    chipClass:
      'border-amber-200 bg-amber-50 text-amber-800 dark:border-amber-700/60 dark:bg-amber-950/35 dark:text-amber-200',
    rowClass:
      'bg-amber-500/[0.05] hover:bg-amber-500/[0.10] dark:bg-amber-400/[0.06] dark:hover:bg-amber-400/[0.11]',
    iconClass: 'bg-amber-500/15 text-amber-600 dark:bg-amber-400/15 dark:text-amber-300',
  },
};

export function canonicalEntityType(value: string | null | undefined): string | undefined {
  return value ? ENTITY_LOOKUP[value.trim().toLowerCase()] : undefined;
}

export function getEntityTypeEmoji(value: string): string {
  const key = canonicalEntityType(value);
  return (key && entityTypeEmojis[key]) || '🔹';
}

export function getEntityTypeStyle(value: string | null | undefined): EntityTypeStyle {
  const key = canonicalEntityType(value);
  const style = (key && entityTypeStyles[key]) || defaultEntityTypeStyle;
  return { ...style, label: (key && BIOLINK_ENTITIES[key]?.label) || style.label };
}
