import 'server-only';

import type { NomenclatureExplanation } from '@/shared/types';

const MOLECULENS_API_URL = process.env.MOLECULENS_API_URL || 'https://api.moleculens.com';
const DEFAULT_TIMEOUT_MS = 10_000;

export class NomenclatureUpstreamError extends Error {
  constructor(
    message: string,
    readonly status: number,
    readonly body: string
  ) {
    super(message);
    this.name = 'NomenclatureUpstreamError';
  }
}

export async function fetchNomenclatureExplanation(input: {
  smiles: string;
  molblock?: string | null;
  verifyOpsin?: boolean;
  timeoutMs?: number;
}): Promise<NomenclatureExplanation> {
  const controller = new AbortController();
  const timeout = setTimeout(() => controller.abort(), input.timeoutMs ?? DEFAULT_TIMEOUT_MS);

  try {
    const response = await fetch(`${MOLECULENS_API_URL}/v1/nomenclature/explain`, {
      method: 'POST',
      headers: { 'content-type': 'application/json', accept: 'application/json' },
      body: JSON.stringify({
        smiles: input.smiles,
        ...(input.molblock ? { molblock: input.molblock } : {}),
        verifyOpsin: input.verifyOpsin ?? true,
      }),
      signal: controller.signal,
      cache: 'no-store',
    });
    const text = await response.text();
    if (!response.ok) {
      throw new NomenclatureUpstreamError(
        `Nomenclature service returned ${response.status}`,
        response.status,
        text.slice(0, 1000)
      );
    }
    if (!text.trim()) {
      throw new NomenclatureUpstreamError('Nomenclature service returned no data', 502, '');
    }
    return JSON.parse(text) as NomenclatureExplanation;
  } catch (error) {
    if (error instanceof Error && error.name === 'AbortError') {
      throw new NomenclatureUpstreamError('Nomenclature service timeout', 504, '');
    }
    throw error;
  } finally {
    clearTimeout(timeout);
  }
}
