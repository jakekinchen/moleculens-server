/**
 * Electrostatics route schema and response transformers.
 */
import { z } from 'zod';

/**
 * Request schema for electrostatics computation.
 */
export const electrostaticsRequestSchema = z.object({
  sdfContent: z.string().min(1, 'sdfContent is required'),
  clientCacheKey: z.string().trim().min(1).optional(),
  method: z.string().optional().default('gasteiger'),
  charge: z.number().int().optional().default(0),
  multiplicity: z.number().int().positive().optional().default(1),
  inchiKey: z.string().optional(),
  presentationId: z.string().uuid().optional(),
  tourHash: z.string().trim().min(1).max(256).optional(),
  no_reorient: z.boolean().optional().default(true),
  no_com: z.boolean().optional().default(true),
  surface: z
    .object({
      gridSpacing: z.number().positive().optional().default(0.25),
      probeRadius: z.number().positive().optional().default(0.9),
      padding: z.number().positive().optional().default(1.2),
    })
    .optional(),
  potential: z
    .object({
      units: z.string().optional().default('kcal/mol/e'),
      softeningEpsilon: z.number().positive().optional().default(0.3),
      clampPercentiles: z.tuple([z.number(), z.number()]).optional().default([5.0, 95.0]),
    })
    .optional(),
});

export type ElectrostaticsRequest = z.infer<typeof electrostaticsRequestSchema>;

/**
 * Build upstream request body from validated payload.
 * Strips fields used only for access control or client-side operations:
 * - presentationId, tourHash: access control only
 * - inchiKey: molecular identity - NOT sent to avoid pose-mismatched cache hits
 */
export function buildUpstreamRequest(payload: ElectrostaticsRequest): Record<string, unknown> {
  const {
    presentationId: _presentationId,
    tourHash: _tourHash,
    inchiKey: _inchiKey,
    clientCacheKey,
    ...forwardPayload
  } = payload;

  const requestBody: Record<string, unknown> = {
    ...forwardPayload,
  };

  if (clientCacheKey) requestBody.client_cache_key = clientCacheKey;

  return requestBody;
}

/**
 * Extract prefetch response (minimal subset).
 */
export function extractPrefetchResponse(data: unknown): unknown {
  if (!data || typeof data !== 'object') return data;

  const obj = data as Record<string, unknown>;
  return {
    cached: obj.cached ?? false,
    status: obj.status,
    jobId: obj.jobId,
    cacheKey: obj.cacheKey,
    error: obj.error,
  };
}
