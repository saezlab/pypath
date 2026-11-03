import type { Declared, WireFilters } from '../api/contracts';

export interface CvTermReference {
  id: string;
  name: string;
}

export interface SearchFilters extends Partial<Declared<WireFilters>> {
  member_a_id?: string | number;
  member_b_id?: string | number;
  is_directed?: boolean | null;
  signs?: Array<-1 | 0 | 1>;

  license_cv?: string[];
  update_category_cv?: string[];
  content_category_cv_terms?: string[];
  ontology_terms?: string[];
  ontology_ids?: string[];

  parent_entity_ids?: Array<string | number>;
  member_entity_ids?: Array<string | number>;
  parent_entity_types?: string[];
  member_entity_types?: string[];
  association_annotation_terms?: string[];

  include_associated_entities?: boolean;
  include_members_participants?: boolean;
  selection_scope_mode?: 'union' | 'intersection';
  scope_endpoint_mode?: 'any' | 'both';
  scope_annotation_ids?: string[];
}

export interface SearchParams {
  query: string;
  filters: SearchFilters;
  limit: number;
  offset: number;
}

export interface SourceFunctionRecord {
  function: string;
  records: number;
}

export interface SearchSource {
  __doc_id?: string;
  content_hash?: string;
  source_ref: string;
  source: string;
  source_name: string;
  source_accession: string;
  license_cv: string;
  update_category_cv: string;
  resource_url?: string;
  resource_description?: string;
  pubmed: string[];
  finished_at: string;
  function_names: string[];
  content_category_cv_terms?: string[];
  function_records: SourceFunctionRecord[];
  total_records: number;
  function_records_json?: string;
  [key: string]: unknown;
}
