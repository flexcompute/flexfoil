import { create } from 'zustand';

export type CaseLogKind = 'single-alpha' | 'single-cl' | 'polar';
export type CaseLogStatus = 'running' | 'success' | 'warning' | 'error';
export type CaseLogLevel = 'info' | 'success' | 'warning' | 'error';

export interface CaseLogEvent {
  id: string;
  timestamp: string;
  elapsedMs: number;
  level: CaseLogLevel;
  message: string;
  details?: Record<string, unknown>;
}

export interface LoggedCase {
  id: string;
  kind: CaseLogKind;
  title: string;
  status: CaseLogStatus;
  startedAt: string;
  startedAtMs: number;
  finishedAt: string | null;
  durationMs: number | null;
  summary: string | null;
  metadata: Record<string, unknown>;
  events: CaseLogEvent[];
}

interface StartCaseInput {
  kind: CaseLogKind;
  title: string;
  metadata?: Record<string, unknown>;
}

interface AppendEventInput {
  level?: CaseLogLevel;
  message: string;
  details?: Record<string, unknown>;
}

interface FinishCaseInput {
  status: Exclude<CaseLogStatus, 'running'>;
  summary?: string;
  metadata?: Record<string, unknown>;
}

interface CaseLogStoreState {
  cases: LoggedCase[];
  selectedCaseId: string | null;
  startCase: (input: StartCaseInput) => string;
  appendEvent: (caseId: string, event: AppendEventInput) => void;
  finishCase: (caseId: string, input: FinishCaseInput) => void;
  selectCase: (caseId: string | null) => void;
  clearCases: () => void;
}

const MAX_CASES = 100;
const MAX_EVENTS_PER_CASE = 400;

function makeId(prefix: string): string {
  return `${prefix}_${Date.now().toString(36)}_${Math.random().toString(36).slice(2, 8)}`;
}

function trimCases(cases: LoggedCase[]): LoggedCase[] {
  return cases.slice(0, MAX_CASES);
}

function trimEvents(events: CaseLogEvent[]): CaseLogEvent[] {
  if (events.length <= MAX_EVENTS_PER_CASE) {
    return events;
  }
  return events.slice(events.length - MAX_EVENTS_PER_CASE);
}

export const useCaseLogStore = create<CaseLogStoreState>()((set, get) => ({
  cases: [],
  selectedCaseId: null,

  startCase: ({ kind, title, metadata = {} }) => {
    const startedAtMs = Date.now();
    const entry: LoggedCase = {
      id: makeId('case'),
      kind,
      title,
      status: 'running',
      startedAt: new Date(startedAtMs).toISOString(),
      startedAtMs,
      finishedAt: null,
      durationMs: null,
      summary: null,
      metadata,
      events: [],
    };

    set((state) => ({
      cases: trimCases([entry, ...state.cases]),
      selectedCaseId: entry.id,
    }));

    return entry.id;
  },

  appendEvent: (caseId, event) => {
    const now = Date.now();
    set((state) => ({
      cases: state.cases.map((item) => {
        if (item.id !== caseId) return item;

        const nextEvent: CaseLogEvent = {
          id: makeId('evt'),
          timestamp: new Date(now).toISOString(),
          elapsedMs: Math.max(0, now - item.startedAtMs),
          level: event.level ?? 'info',
          message: event.message,
          details: event.details,
        };

        return {
          ...item,
          events: trimEvents([...item.events, nextEvent]),
        };
      }),
    }));
  },

  finishCase: (caseId, input) => {
    const now = Date.now();
    set((state) => ({
      cases: state.cases.map((item) => {
        if (item.id !== caseId) return item;
        return {
          ...item,
          status: input.status,
          finishedAt: new Date(now).toISOString(),
          durationMs: Math.max(0, now - item.startedAtMs),
          summary: input.summary ?? item.summary,
          metadata: input.metadata ? { ...item.metadata, ...input.metadata } : item.metadata,
        };
      }),
      selectedCaseId: get().selectedCaseId ?? caseId,
    }));
  },

  selectCase: (caseId) => set({ selectedCaseId: caseId }),

  clearCases: () => set({ cases: [], selectedCaseId: null }),
}));
