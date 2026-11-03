export type MolecularFormRecord = {
  protein_entity_key?: string | null;
  transcript_entity_key?: string | null;
  isoform_identifier?: { ns: string; id: string } | null;
  sequence_identifiers?: Array<{ ns: string; id: string }> | null;
  modifications?: Array<Record<string, unknown>> | null;
  variants?: Array<Record<string, unknown>> | null;
  regions?: Array<Record<string, unknown>> | null;
};
