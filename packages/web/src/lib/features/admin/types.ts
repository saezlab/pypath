import type { Declared, WireAdminJob, WireAdminStage } from '../../api/contracts';

export type LogInfo = {
  exists: boolean;
  size_bytes?: number;
  mtime?: number | null;
};

export type FileInfo = {
  name: string;
  exists: boolean;
  rows: number;
  size_bytes: number;
  mtime: number | null;
  log?: LogInfo;
};

export type ResolutionCounts = {
  input: number;
  resolved: number;
  unresolved: number;
  not_applicable: number;
};

export type ResolutionStats = {
  scope?: string;
  input_entities: number;
  resolved_entities: number;
  /** Chemicals identified by a valid InChIKey that the reference does not hold yet. */
  structure_entities?: number;
  unresolved_entities: number;
  not_applicable_entities: number;
  lookup_entities?: number;
  by_entity_type?: Record<string, ResolutionCounts>;
};

export type Stage = Omit<Declared<WireAdminStage>, 'current' | 'total'> & {
  current: number | null;
  total: number | null;
};
export type Job = Omit<Declared<WireAdminJob>, 'params' | 'stages' | 'logs' | 'error'> & {
  params: Record<string, unknown>;
  stages: Stage[];
  logs: string[];
  error: string | null;
};

export type BuiltResource = {
  resource: string;
  version: string;
  entities: FileInfo;
  relations: FileInfo;
  payloads: FileInfo;
  resolution_stats?: ResolutionStats | null;
  log?: LogInfo;
};

export type AdminStatus = {
  build_available: boolean;
  hubs_available?: boolean;
  pipeline_available?: boolean;
  resolver: {
    ready: boolean;
    hubs: FileInfo[];
  };
  resources: { built: BuiltResource[]; logs?: Record<string, LogInfo> };
  job: Job | null;
  jobs?: Job[];
};

export type AvailableSource = { source: string; datasets: string[] };
