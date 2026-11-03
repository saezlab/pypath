import type { EntityRelation, Entity } from '$lib/domain/entity-types';

export type InteractionListRow = {
  relation: EntityRelation;
  subjectEntity: Entity;
  objectEntity: Entity;
};

export type InteractionDetailsData = InteractionListRow & {
  evidence: Array<{
    relationPk: string;
    source: string;
    subjectMolecularForm?: import('$lib/types/molecular').MolecularFormRecord | null;
    objectMolecularForm?: import('$lib/types/molecular').MolecularFormRecord | null;
    subjectAttributes: unknown;
    recordAttributes: unknown;
    evidence: unknown;
    objectAttributes: unknown;
  }>;
};

export type ParsedAnnotation = {
  term: string;
  termId?: string;
  value?: string;
  unit?: string;
  unitId?: string;
  sourceField?: string;
};

export type EvidenceGroup = {
  key: string;
  source: string;
  evidenceCount: number;
  pubmedIds: string[];
  subjectAnnotations: ParsedAnnotation[];
  relationAnnotations: ParsedAnnotation[];
  objectAnnotations: ParsedAnnotation[];
};
