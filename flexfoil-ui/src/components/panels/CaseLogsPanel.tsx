import { useMemo } from 'react';
import {
  useCaseLogStore,
  type CaseLogLevel,
  type CaseLogStatus,
  type LoggedCase,
} from '../../stores/caseLogStore';

function formatElapsed(elapsedMs: number): string {
  if (elapsedMs < 1000) return `${elapsedMs} ms`;
  return `${(elapsedMs / 1000).toFixed(2)} s`;
}

function formatDateTime(value: string | null): string {
  if (!value) return 'In progress';
  return new Date(value).toLocaleTimeString();
}

function statusColor(status: CaseLogStatus): string {
  switch (status) {
    case 'success':
      return 'var(--accent-primary)';
    case 'warning':
      return 'var(--accent-warning)';
    case 'error':
      return 'var(--accent-danger)';
    default:
      return 'var(--accent-secondary)';
  }
}

function levelColor(level: CaseLogLevel): string {
  switch (level) {
    case 'success':
      return 'var(--accent-primary)';
    case 'warning':
      return 'var(--accent-warning)';
    case 'error':
      return 'var(--accent-danger)';
    default:
      return 'var(--text-secondary)';
  }
}

function kindLabel(kind: LoggedCase['kind']): string {
  switch (kind) {
    case 'single-alpha':
      return 'Single α';
    case 'single-cl':
      return 'Single CL';
    case 'polar':
      return 'Polar';
  }
}

function renderValue(value: unknown): string {
  if (typeof value === 'number') {
    if (!Number.isFinite(value)) return String(value);
    if (Math.abs(value) >= 1000 || (Math.abs(value) > 0 && Math.abs(value) < 1e-3)) {
      return value.toExponential(3);
    }
    return Number.isInteger(value) ? String(value) : value.toFixed(4).replace(/\.?0+$/, '');
  }
  if (typeof value === 'boolean') return value ? 'true' : 'false';
  if (value == null) return 'null';
  if (typeof value === 'string') return value;
  return JSON.stringify(value);
}

function DetailGrid({ data }: { data: Record<string, unknown> }) {
  const entries = Object.entries(data);
  if (entries.length === 0) {
    return (
      <div style={{ fontSize: '11px', color: 'var(--text-muted)' }}>
        No metadata recorded.
      </div>
    );
  }

  return (
    <div
      style={{
        display: 'grid',
        gridTemplateColumns: 'max-content 1fr',
        gap: '6px 10px',
        fontSize: '11px',
        alignItems: 'start',
      }}
    >
      {entries.map(([key, value]) => (
        <FragmentRow key={key} label={key} value={renderValue(value)} />
      ))}
    </div>
  );
}

function FragmentRow({ label, value }: { label: string; value: string }) {
  return (
    <>
      <div style={{ color: 'var(--text-muted)', textTransform: 'capitalize' }}>
        {label.replace(/_/g, ' ')}
      </div>
      <div style={{ fontFamily: 'var(--font-mono)', wordBreak: 'break-word' }}>{value}</div>
    </>
  );
}

export function CaseLogsPanel() {
  const cases = useCaseLogStore((state) => state.cases);
  const selectedCaseId = useCaseLogStore((state) => state.selectedCaseId);
  const selectCase = useCaseLogStore((state) => state.selectCase);
  const clearCases = useCaseLogStore((state) => state.clearCases);

  const selectedCase = useMemo(
    () => cases.find((entry) => entry.id === selectedCaseId) ?? cases[0] ?? null,
    [cases, selectedCaseId],
  );

  return (
    <div className="panel">
      <div className="panel-header">Case Logs</div>
      <div
        className="panel-content"
        style={{ padding: 0, display: 'flex', flexDirection: 'column', overflow: 'hidden' }}
      >
        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            gap: '8px',
            padding: '10px 12px',
            borderBottom: '1px solid var(--border-color)',
            flexShrink: 0,
          }}
        >
          <div style={{ fontSize: '11px', color: 'var(--text-muted)' }}>
            {cases.length} captured case{cases.length === 1 ? '' : 's'}
          </div>
          <div style={{ flex: 1 }} />
          <button
            onClick={clearCases}
            disabled={cases.length === 0}
            style={{ padding: '4px 8px', fontSize: '11px' }}
          >
            Clear
          </button>
        </div>

        {cases.length === 0 ? (
          <div
            style={{
              flex: 1,
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              padding: '16px',
              fontSize: '12px',
              color: 'var(--text-muted)',
              textAlign: 'center',
            }}
          >
            Run a single-point solve or polar sweep to capture detailed case logs.
          </div>
        ) : (
          <div style={{ flex: 1, minHeight: 0, display: 'flex' }}>
            <div
              style={{
                width: '34%',
                minWidth: '220px',
                borderRight: '1px solid var(--border-color)',
                overflowY: 'auto',
              }}
            >
              {cases.map((entry) => {
                const active = entry.id === selectedCase?.id;
                return (
                  <button
                    key={entry.id}
                    onClick={() => selectCase(entry.id)}
                    style={{
                      width: '100%',
                      textAlign: 'left',
                      border: 'none',
                      borderBottom: '1px solid var(--border-color)',
                      borderRadius: 0,
                      padding: '10px 12px',
                      background: active ? 'var(--bg-hover)' : 'transparent',
                    }}
                  >
                    <div
                      style={{
                        display: 'flex',
                        alignItems: 'center',
                        gap: '8px',
                        marginBottom: '6px',
                      }}
                    >
                      <span
                        style={{
                          fontSize: '10px',
                          fontWeight: 600,
                          color: statusColor(entry.status),
                          textTransform: 'uppercase',
                          letterSpacing: '0.4px',
                        }}
                      >
                        {entry.status}
                      </span>
                      <span
                        style={{
                          fontSize: '10px',
                          color: 'var(--text-muted)',
                          padding: '2px 6px',
                          border: '1px solid var(--border-color)',
                          borderRadius: '999px',
                        }}
                      >
                        {kindLabel(entry.kind)}
                      </span>
                    </div>
                    <div style={{ fontSize: '12px', fontWeight: 600, marginBottom: '4px' }}>
                      {entry.title}
                    </div>
                    <div style={{ fontSize: '10px', color: 'var(--text-muted)' }}>
                      {formatDateTime(entry.startedAt)}
                      {entry.durationMs != null ? ` • ${formatElapsed(entry.durationMs)}` : ''}
                    </div>
                    {entry.summary && (
                      <div
                        style={{
                          marginTop: '6px',
                          fontSize: '10px',
                          color: 'var(--text-secondary)',
                          lineHeight: 1.4,
                        }}
                      >
                        {entry.summary}
                      </div>
                    )}
                  </button>
                );
              })}
            </div>

            <div style={{ flex: 1, minWidth: 0, display: 'flex', flexDirection: 'column' }}>
              {selectedCase && (
                <>
                  <div
                    style={{
                      padding: '12px',
                      borderBottom: '1px solid var(--border-color)',
                      flexShrink: 0,
                    }}
                  >
                    <div
                      style={{
                        display: 'flex',
                        alignItems: 'center',
                        gap: '8px',
                        marginBottom: '8px',
                        flexWrap: 'wrap',
                      }}
                    >
                      <div style={{ fontSize: '14px', fontWeight: 700 }}>{selectedCase.title}</div>
                      <span
                        style={{
                          fontSize: '10px',
                          color: statusColor(selectedCase.status),
                          textTransform: 'uppercase',
                          letterSpacing: '0.5px',
                        }}
                      >
                        {selectedCase.status}
                      </span>
                    </div>

                    <div style={{ fontSize: '11px', color: 'var(--text-muted)', marginBottom: '10px' }}>
                      Started {formatDateTime(selectedCase.startedAt)}
                      {selectedCase.finishedAt ? ` • Finished ${formatDateTime(selectedCase.finishedAt)}` : ''}
                      {selectedCase.durationMs != null ? ` • Duration ${formatElapsed(selectedCase.durationMs)}` : ''}
                    </div>

                    {selectedCase.summary && (
                      <div
                        style={{
                          marginBottom: '12px',
                          padding: '8px 10px',
                          background: 'var(--bg-tertiary)',
                          borderRadius: '6px',
                          fontSize: '11px',
                          lineHeight: 1.5,
                        }}
                      >
                        {selectedCase.summary}
                      </div>
                    )}

                    <DetailGrid data={selectedCase.metadata} />
                  </div>

                  <div style={{ flex: 1, minHeight: 0, overflowY: 'auto', padding: '12px' }}>
                    <div
                      style={{
                        fontSize: '11px',
                        textTransform: 'uppercase',
                        letterSpacing: '0.5px',
                        color: 'var(--text-muted)',
                        marginBottom: '10px',
                      }}
                    >
                      Event Timeline
                    </div>

                    {selectedCase.events.length === 0 ? (
                      <div style={{ fontSize: '11px', color: 'var(--text-muted)' }}>
                        No events recorded yet.
                      </div>
                    ) : (
                      <div style={{ display: 'flex', flexDirection: 'column', gap: '8px' }}>
                        {selectedCase.events.map((event) => (
                          <div
                            key={event.id}
                            style={{
                              padding: '8px 10px',
                              border: '1px solid var(--border-color)',
                              borderLeft: `3px solid ${levelColor(event.level)}`,
                              borderRadius: '6px',
                              background: 'var(--bg-tertiary)',
                            }}
                          >
                            <div
                              style={{
                                display: 'flex',
                                alignItems: 'center',
                                gap: '8px',
                                flexWrap: 'wrap',
                                marginBottom: event.details ? '6px' : 0,
                              }}
                            >
                              <span
                                style={{
                                  fontSize: '10px',
                                  color: levelColor(event.level),
                                  textTransform: 'uppercase',
                                  letterSpacing: '0.4px',
                                }}
                              >
                                {event.level}
                              </span>
                              <span style={{ fontSize: '10px', color: 'var(--text-muted)' }}>
                                +{formatElapsed(event.elapsedMs)}
                              </span>
                              <span style={{ fontSize: '11px', color: 'var(--text-primary)' }}>
                                {event.message}
                              </span>
                            </div>

                            {event.details && (
                              <pre
                                style={{
                                  margin: 0,
                                  fontSize: '10px',
                                  whiteSpace: 'pre-wrap',
                                  wordBreak: 'break-word',
                                  fontFamily: 'var(--font-mono)',
                                  color: 'var(--text-secondary)',
                                }}
                              >
                                {JSON.stringify(event.details, null, 2)}
                              </pre>
                            )}
                          </div>
                        ))}
                      </div>
                    )}
                  </div>
                </>
              )}
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
