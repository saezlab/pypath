import {
  Activity,
  ArrowRightLeft,
  CircleDot,
  Dna,
  FileText,
  FlaskConical,
  HelpCircle,
  Leaf,
  Route,
  Shapes,
  Tag,
  Users,
  Waves,
} from '@lucide/svelte';
import { canonicalEntityType } from './entity-types';

type Icon = typeof CircleDot;

// Keyed by canonical Biolink class name.
const ICONS: Record<string, Icon> = {
  gene: Dna,
  protein: CircleDot,
  polypeptide: CircleDot,
  protein_family: Users,
  gene_family: Users,
  chemical_entity: FlaskConical,
  small_molecule: FlaskConical,
  molecular_entity: FlaskConical,
  food: FlaskConical,
  macromolecular_complex: Shapes,
  physical_entity: Shapes,
  pathway: Route,
  molecular_activity: ArrowRightLeft,
  biological_process: Activity,
  ontology_class: Tag,
  information_content_entity: FileText,
  microrna: Waves,
  rna_product: Waves,
  organism_taxon: Leaf,
};

// Display labels and legacy names that are not Biolink class names.
const ALIASES: Record<string, string> = {
  reaction: 'molecular_activity',
  chemical: 'chemical_entity',
  compound: 'chemical_entity',
  metabolite: 'chemical_entity',
  drug: 'chemical_entity',
  lipid: 'chemical_entity',
  complex: 'macromolecular_complex',
  mirna: 'microrna',
  cvterm: 'ontology_class',
  ontologyterm: 'ontology_class',
  enzymeclass: 'ontology_class',
  keyword: 'ontology_class',
  chemicalclass: 'ontology_class',
  figure: 'information_content_entity',
  informationrecord: 'information_content_entity',
};

/** Lucide icon for an entity type, given its Biolink class name or display label. */
export function getEntityTypeIcon(value: string | null | undefined): Icon {
  if (!value) return HelpCircle;
  const bare = value.replace(/^biolink:/i, '');
  const typeName = bare.includes(':') ? bare.split(':')[0] : bare;
  const key =
    canonicalEntityType(typeName) ?? ALIASES[typeName.toLowerCase().replace(/[\s_-]/g, '')];
  return (key && ICONS[key]) || HelpCircle;
}
