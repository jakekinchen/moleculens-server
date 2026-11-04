import type { ElectrostaticsComputeRequest } from '@/shared/lib/electrostatics/types';
import type { OrbitalComputeRequest, OrbitalType } from '@/shared/lib/orbitals/types';

const OVERLAY_GEOMETRY_CACHE_KEY_VERSION = 'overlay-geometry-v1';
const OVERLAY_CLIENT_CACHE_KEY_VERSION = 'overlay-v3';

function fallbackHash(input: string): string {
  let hash = 0;
  for (let i = 0; i < input.length; i++) {
    const char = input.charCodeAt(i);
    hash = ((hash << 5) - hash + char) | 0;
  }
  return Math.abs(hash).toString(16);
}

function sortOrbitals(orbitals: OrbitalType[]): OrbitalType[] {
  return [...new Set(orbitals)].sort();
}

function normalizeForHash(value: unknown): unknown {
  if (Array.isArray(value)) {
    return value.map(item => normalizeForHash(item));
  }

  if (value && typeof value === 'object') {
    return Object.fromEntries(
      Object.entries(value as Record<string, unknown>)
        .filter(([, entryValue]) => entryValue !== undefined)
        .sort(([left], [right]) => left.localeCompare(right))
        .map(([key, entryValue]) => [key, normalizeForHash(entryValue)])
    );
  }

  return value;
}

async function hashPayload(payload: string): Promise<string> {
  if (typeof crypto !== 'undefined' && crypto.subtle) {
    const encoder = new TextEncoder();
    const data = encoder.encode(payload);
    const hashBuffer = await crypto.subtle.digest('SHA-256', data);
    const hashArray = Array.from(new Uint8Array(hashBuffer));
    return hashArray.map(byte => byte.toString(16).padStart(2, '0')).join('');
  }

  return fallbackHash(payload);
}

export async function computeOverlayGeometryCacheKey(sdf: string): Promise<string> {
  const payload = `${OVERLAY_GEOMETRY_CACHE_KEY_VERSION}:${sdf}`;
  return hashPayload(payload);
}

export async function computeOrbitalClientCacheKey(
  args: Pick<
    OrbitalComputeRequest,
    | 'sdf'
    | 'orbitals'
    | 'charge'
    | 'multiplicity'
    | 'method'
    | 'basis'
    | 'gridSpacing'
    | 'isovalue'
    | 'densityIsovalue'
  >
): Promise<string> {
  const payload = JSON.stringify(
    normalizeForHash({
      version: OVERLAY_CLIENT_CACHE_KEY_VERSION,
      kind: 'orbitals',
      sdf: args.sdf,
      orbitals: sortOrbitals(args.orbitals),
      charge: args.charge ?? 0,
      multiplicity: args.multiplicity ?? 1,
      method: args.method ?? 'scf',
      basis: args.basis ?? 'sto-3g',
      gridSpacing: args.gridSpacing ?? 0.25,
      isovalue: args.isovalue ?? 0.02,
      densityIsovalue: args.densityIsovalue ?? 0.01,
    })
  );

  return hashPayload(payload);
}

export async function computeElectrostaticsClientCacheKey(
  args: Pick<
    ElectrostaticsComputeRequest,
    | 'sdfContent'
    | 'method'
    | 'charge'
    | 'multiplicity'
    | 'no_reorient'
    | 'no_com'
    | 'surface'
    | 'potential'
  >
): Promise<string> {
  const payload = JSON.stringify(
    normalizeForHash({
      version: OVERLAY_CLIENT_CACHE_KEY_VERSION,
      kind: 'electrostatics',
      sdfContent: args.sdfContent,
      method: args.method ?? 'gasteiger',
      charge: args.charge ?? 0,
      multiplicity: args.multiplicity ?? 1,
      no_reorient: args.no_reorient ?? true,
      no_com: args.no_com ?? true,
      surface: {
        gridSpacing: args.surface?.gridSpacing ?? 0.25,
        probeRadius: args.surface?.probeRadius ?? 0.9,
        padding: args.surface?.padding ?? 1.2,
      },
      potential: {
        units: args.potential?.units ?? 'kcal/mol/e',
        softeningEpsilon: args.potential?.softeningEpsilon ?? 0.3,
        clampPercentiles: args.potential?.clampPercentiles ?? [5.0, 95.0],
      },
    })
  );

  return hashPayload(payload);
}
