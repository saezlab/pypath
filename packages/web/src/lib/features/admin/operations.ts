import type { WireAdminJob, WireAdminJobs } from '../../api/contracts';
import { adminFetch } from './auth';
import type { AdminStatus, AvailableSource, Job } from './types';

export class AdminRequestError extends Error {
  constructor(
    public status: number,
    message: string,
  ) {
    super(message);
  }
}
async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await adminFetch(`/app-api/admin/${path}`, init);
  if (!response.ok) {
    const body = (await response.json().catch(() => ({}))) as { detail?: string; error?: string };
    throw new AdminRequestError(
      response.status,
      String(body.detail || body.error || response.statusText),
    );
  }
  return response.json() as Promise<T>;
}
export function jobFromWire(job: WireAdminJob): Job {
  return {
    ...job,
    params: job.params ?? {},
    logs: job.logs ?? [],
    error: job.error ?? null,
    stages: (job.stages ?? []).map((stage) => ({
      ...stage,
      current: stage.current ?? null,
      total: stage.total ?? null,
    })),
  };
}
export async function fetchAdminStatus(signal?: AbortSignal): Promise<AdminStatus> {
  const data = await request<Omit<AdminStatus, 'job' | 'jobs'> & WireAdminJobs>('status', {
    signal,
  });
  return {
    ...data,
    job: data.job ? jobFromWire(data.job) : null,
    jobs: (data.jobs ?? []).map(jobFromWire),
  };
}
export const fetchAdminSources = async () =>
  (await request<{ sources: AvailableSource[] }>('sources')).sources;
export const submitAdminJob = (body: Record<string, unknown>) =>
  request('jobs', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });
export const cancelAdminJob = (id: string) =>
  request(`jobs/${encodeURIComponent(id)}/cancel`, { method: 'POST' });
export const removeResourceVersion = (resource: string, version: string) =>
  request(`resources/${encodeURIComponent(resource)}/versions/${encodeURIComponent(version)}`, {
    method: 'DELETE',
  });
export async function fetchAdminLog(
  kind: string,
  name: string,
  version: string,
  signal?: AbortSignal,
) {
  const query = version ? `?version=${encodeURIComponent(version)}` : '';
  return (
    await request<{ text: string }>(
      `logs/${encodeURIComponent(kind)}/${encodeURIComponent(name)}${query}`,
      { signal },
    )
  ).text;
}
