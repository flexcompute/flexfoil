/**
 * Zustand store tracking solver jobs for the global status bar.
 *
 * Call `dispatch` before starting a solver call.  It returns an `id`
 * and an `AbortSignal`.  Check `signal.aborted` before/between
 * expensive WASM calls.  Call `complete` or `cancel` when done.
 */

import { create } from 'zustand';

export interface SolverJob {
  id: string;
  label: string;
  status: 'running' | 'done' | 'cancelled' | 'error';
  startedAt: number;
  finishedAt?: number;
  progress?: string;
  error?: string;
  _abort?: AbortController;
}

interface SolverJobStore {
  jobs: SolverJob[];

  /** Start a new job.  Returns an id + AbortSignal for cancellation. */
  dispatch: (label: string) => { id: string; signal: AbortSignal };

  /** Update progress text on a running job. */
  update: (id: string, progress: string) => void;

  /** Mark a job as finished (success or error). */
  complete: (id: string, error?: string) => void;

  /** Cancel a running job (sets signal.aborted = true). */
  cancel: (id: string) => void;

  /** Remove all finished/cancelled jobs from the list. */
  clearHistory: () => void;
}

let nextId = 1;

export const useSolverJobStore = create<SolverJobStore>((set, get) => ({
  jobs: [],

  dispatch: (label) => {
    const id = `job-${nextId++}`;
    const abort = new AbortController();
    const job: SolverJob = {
      id,
      label,
      status: 'running',
      startedAt: Date.now(),
      _abort: abort,
    };
    set((s) => ({ jobs: [...s.jobs, job] }));
    return { id, signal: abort.signal };
  },

  update: (id, progress) =>
    set((s) => ({
      jobs: s.jobs.map((j) => (j.id === id ? { ...j, progress } : j)),
    })),

  complete: (id, error) =>
    set((s) => ({
      jobs: s.jobs.map((j) =>
        j.id === id
          ? {
              ...j,
              status: error ? 'error' as const : 'done' as const,
              finishedAt: Date.now(),
              error,
            }
          : j,
      ),
    })),

  cancel: (id) => {
    const job = get().jobs.find((j) => j.id === id);
    if (job?._abort) job._abort.abort();
    set((s) => ({
      jobs: s.jobs.map((j) =>
        j.id === id
          ? { ...j, status: 'cancelled' as const, finishedAt: Date.now() }
          : j,
      ),
    }));
  },

  clearHistory: () =>
    set((s) => ({
      jobs: s.jobs.filter((j) => j.status === 'running'),
    })),
}));

/** Convenience: get the currently running job (if any). */
export function useActiveJob(): SolverJob | null {
  return useSolverJobStore((s) => s.jobs.find((j) => j.status === 'running') ?? null);
}
