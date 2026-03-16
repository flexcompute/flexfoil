import { beforeEach, describe, expect, it } from 'vitest';
import { useCaseLogStore } from './caseLogStore';

describe('caseLogStore', () => {
  beforeEach(() => {
    useCaseLogStore.getState().clearCases();
  });

  it('starts a case and selects it automatically', () => {
    const caseId = useCaseLogStore.getState().startCase({
      kind: 'single-alpha',
      title: 'Run to α 4.00°',
      metadata: { airfoil: 'NACA 0012' },
    });

    const state = useCaseLogStore.getState();
    expect(state.cases).toHaveLength(1);
    expect(state.selectedCaseId).toBe(caseId);
    expect(state.cases[0].metadata.airfoil).toBe('NACA 0012');
    expect(state.cases[0].status).toBe('running');
  });

  it('appends events and finishes a case with summary data', () => {
    const caseId = useCaseLogStore.getState().startCase({
      kind: 'polar',
      title: 'Alpha polar -5.0° to 15.0°',
    });

    useCaseLogStore.getState().appendEvent(caseId, {
      message: 'Starting polar sweep',
      details: { alpha_start: -5, alpha_end: 15, alpha_step: 1 },
    });
    useCaseLogStore.getState().finishCase(caseId, {
      status: 'success',
      summary: 'Polar sweep completed with 21 valid points.',
      metadata: { cache_hits: 3, computed_points: 18 },
    });

    const entry = useCaseLogStore.getState().cases[0];
    expect(entry.events).toHaveLength(1);
    expect(entry.events[0].message).toBe('Starting polar sweep');
    expect(entry.summary).toContain('21 valid points');
    expect(entry.status).toBe('success');
    expect(entry.metadata.cache_hits).toBe(3);
    expect(entry.finishedAt).not.toBeNull();
  });

  it('clears all captured cases', () => {
    useCaseLogStore.getState().startCase({
      kind: 'single-cl',
      title: 'Run to CL 0.500',
    });

    useCaseLogStore.getState().clearCases();

    const state = useCaseLogStore.getState();
    expect(state.cases).toHaveLength(0);
    expect(state.selectedCaseId).toBeNull();
  });
});
