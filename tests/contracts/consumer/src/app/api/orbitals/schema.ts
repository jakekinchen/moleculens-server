/**
 * Orbitals route schema and response transformers.
 */
import { z } from 'zod';

import { ORBITAL_TYPES } from '@/shared/lib/orbitals/types';

/**
 * Request schema for orbital computation.
 */
export const orbitalsRequestSchema = z.object({
  sdf: z.string().min(1, 'SDF is required'),
  smiles: z.string().optional(),
  inchiKey: z.string().optional(),
  clientCacheKey: z.string().optional(),
  presentationId: z.string().uuid().optional(),
  tourHash: z.string().trim().min(1).max(256).optional(),
  charge: z.number().int().optional().default(0),
  multiplicity: z.number().int().positive().optional().default(1),
  method: z.string().optional().default('scf'),
  basis: z.string().optional().default('sto-3g'),
  orbitals: z.array(z.enum(ORBITAL_TYPES)).min(1),
  gridSpacing: z.number().positive().optional().default(0.25),
  isovalue: z.number().positive().optional().default(0.02),
  densityIsovalue: z.number().positive().optional().default(0.01),
});

export type OrbitalsRequest = z.infer<typeof orbitalsRequestSchema>;

/**
 * Build upstream request body from validated payload.
 * Strips loose molecule-identity fields when a recipe-aware client cache key is present:
 * - smiles, inchiKey: can collide across poses/recipes and poison exact-geometry cache entries
 * - presentationId, tourHash: access control only
 */
export function buildUpstreamRequest(payload: OrbitalsRequest): Record<string, unknown> {
  const {
    sdf,
    presentationId: _presentationId,
    tourHash: _tourHash,
    smiles,
    inchiKey,
    clientCacheKey,
    gridSpacing,
    densityIsovalue,
    ...forwardPayload
  } = payload;

  const requestBody: Record<string, unknown> = {
    sdf_content: sdf,
    ...forwardPayload,
    grid_spacing: gridSpacing,
    density_isovalue: densityIsovalue,
  };

  if (clientCacheKey) requestBody.client_cache_key = clientCacheKey;
  if (!clientCacheKey && smiles) requestBody.smiles = smiles;
  if (!clientCacheKey && inchiKey) requestBody.inchi_key = inchiKey;

  return requestBody;
}

interface UpstreamOrbitalResult {
  orbitals?: unknown;
  density?: unknown;
  meta?: {
    method?: string;
    basis?: string;
    grid_spacing_angstrom?: number;
  };
}

interface UpstreamOrbitalError {
  message?: string;
  code?: string;
}

interface UpstreamOrbitalResponse {
  cached?: boolean;
  status?: string;
  cache_key?: string;
  cacheKey?: string;
  job_id?: string;
  jobId?: string;
  created_at?: string;
  createdAt?: string;
  compute_time_ms?: number;
  computeTimeMs?: number;
  progress?: number;
  step?: string;
  result?: UpstreamOrbitalResult;
  error_message?: string;
  errorMessage?: string;
  error?: UpstreamOrbitalError | string;
}

const MISSING_RESULT_ERROR_MESSAGE =
  'Orbital result is unavailable right now. Please try again later.';

/**
 * Transform orbital response from snake_case to camelCase.
 */
export function transformOrbitalsResponse(data: unknown): unknown {
  if (!data || typeof data !== 'object') return data;

  const raw = data as UpstreamOrbitalResponse;
  const result = raw.result
    ? {
        orbitals: raw.result.orbitals,
        density: raw.result.density,
        meta: raw.result.meta
          ? {
              method: raw.result.meta.method,
              basis: raw.result.meta.basis,
              grid_spacing_angstrom: raw.result.meta.grid_spacing_angstrom,
            }
          : undefined,
      }
    : undefined;

  if (raw.status === 'done' && !result) {
    return {
      cached: raw.cached ?? false,
      status: 'error',
      jobId: raw.job_id ?? raw.jobId,
      progress: raw.progress,
      step: raw.step,
      error: {
        message: MISSING_RESULT_ERROR_MESSAGE,
        code: 'missing_result',
      },
    };
  }

  return {
    cached: raw.cached ?? false,
    status: raw.status,
    cacheKey: raw.cache_key ?? raw.cacheKey,
    jobId: raw.job_id ?? raw.jobId,
    createdAt: raw.created_at ?? raw.createdAt,
    computeTimeMs: raw.compute_time_ms ?? raw.computeTimeMs,
    progress: raw.progress,
    step: raw.step,
    result,
    error:
      (raw.error_message ?? raw.errorMessage)
        ? {
            message: raw.error_message ?? raw.errorMessage ?? MISSING_RESULT_ERROR_MESSAGE,
          }
        : typeof raw.error === 'string'
          ? { message: raw.error }
          : raw.error
            ? {
                message: raw.error.message ?? String(raw.error),
                code: raw.error.code,
              }
            : undefined,
  };
}

/**
 * Extract prefetch response (minimal subset).
 */
export function extractPrefetchResponse(data: unknown): unknown {
  if (!data || typeof data !== 'object') return data;

  const obj = data as Record<string, unknown>;
  return {
    cached: obj.cached,
    status: obj.status,
    cacheKey: obj.cacheKey,
    jobId: obj.jobId,
    createdAt: obj.createdAt,
    computeTimeMs: obj.computeTimeMs,
    error: obj.error,
  };
}
