/**
 * Conformers route schema and response transformers.
 */
import { NextResponse } from 'next/server';

/**
 * Request interface for conformer computation.
 */
export interface ConformerRequest {
  molblock2d: string;
  method?: string;
  waitMs?: number;
  seed?: number;
  numConformers?: number;
  optimize?: boolean;
  forceField?: string;
  clientCacheKey?: string;
  smiles?: string;
  inchiKey?: string;
}

/**
 * Build upstream request body from validated payload.
 */
export function buildUpstreamRequest(body: ConformerRequest): Record<string, unknown> {
  return {
    molblock2d: body.molblock2d,
    geomVersion: body.method || 'etkdgv3',
    waitMs: body.waitMs ?? 2000,
    seed: body.seed,
    numConformers: body.numConformers ?? 1,
    optimize: body.optimize ?? true,
    forceField: body.forceField ?? 'mmff',
    clientCacheKey: body.clientCacheKey,
    smiles: body.smiles,
    inchiKey: body.inchiKey,
  };
}

/**
 * Create standardized JSON error response for conformers.
 */
export function jsonError(status: number, payload: Record<string, unknown>) {
  return NextResponse.json(
    {
      error: payload.error || 'error',
      error_code: payload.error_code || 'UNKNOWN',
      message: payload.message || 'An error occurred',
      ...(payload.details ? { details: payload.details } : {}),
    },
    { status }
  );
}

interface UpstreamConformerResponse {
  status?: string;
  cached?: boolean;
  job_id?: string;
  jobId?: string;
  cache_key?: string;
  cacheKey?: string;
  sdf3d?: string;
  sdf_3d?: string;
  progress?: number;
  step?: string;
  error?: unknown;
  errorMessage?: string;
}

/**
 * Transform conformer response from snake_case to camelCase.
 */
export function transformConformersResponse(data: unknown): Record<string, unknown> {
  if (!data || typeof data !== 'object') {
    return { raw: data };
  }

  const obj = data as UpstreamConformerResponse;

  // Normalize error field
  const normalizedError =
    typeof obj.error === 'object' && obj.error !== null
      ? obj.error
      : typeof obj.error === 'string'
        ? { code: 'UPSTREAM_ERROR', message: obj.error }
        : typeof obj.errorMessage === 'string'
          ? { code: 'UPSTREAM_ERROR', message: obj.errorMessage }
          : undefined;

  return {
    status: obj.status,
    cached: obj.cached ?? false,
    jobId: obj.job_id ?? obj.jobId,
    cacheKey: obj.cache_key ?? obj.cacheKey,
    sdf3d: obj.sdf3d ?? obj.sdf_3d,
    progress: obj.progress,
    step: obj.step,
    error: normalizedError,
  };
}
