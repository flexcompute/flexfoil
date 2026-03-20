/**
 * SolverStatusIndicator - Compact status dot+label that sits in the brand footer.
 * Click to expand a popover with recent job history, progress bars, and cancel buttons.
 */

import { useState, useEffect, useRef } from 'react';
import { useSolverJobStore, useActiveJob, type SolverJob } from '../stores/solverJobStore';

function elapsed(ms: number): string {
  if (ms < 1000) return `${ms}ms`;
  const s = ms / 1000;
  if (s < 60) return `${s.toFixed(1)}s`;
  return `${Math.floor(s / 60)}m ${Math.round(s % 60)}s`;
}

function MiniProgressBar({ fraction }: { fraction: number }) {
  return (
    <div style={{
      width: 48, height: 4, borderRadius: 2,
      background: 'var(--bg-tertiary, #333)',
      overflow: 'hidden', flexShrink: 0,
    }}>
      <div style={{
        height: '100%',
        width: `${Math.min(100, fraction * 100)}%`,
        background: 'var(--accent-primary, #3b82f6)',
        borderRadius: 2,
        transition: 'width 0.15s',
      }} />
    </div>
  );
}

function JobRow({ job }: { job: SolverJob }) {
  const cancel = useSolverJobStore((s) => s.cancel);
  const duration = job.finishedAt
    ? elapsed(job.finishedAt - job.startedAt)
    : elapsed(Date.now() - job.startedAt);

  const statusColor =
    job.status === 'running' ? 'var(--accent-warning, #eab308)' :
    job.status === 'done' ? 'var(--accent-success, #22c55e)' :
    job.status === 'error' ? 'var(--accent-error, #ef4444)' :
    'var(--text-muted)';

  return (
    <div style={{
      padding: '5px 0',
      borderBottom: '1px solid var(--border-color)',
      fontSize: 11,
    }}>
      <div style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
        <span style={{
          width: 7, height: 7, borderRadius: '50%', flexShrink: 0,
          background: statusColor,
          animation: job.status === 'running' ? 'solver-pulse 1.2s ease-in-out infinite' : 'none',
        }} />
        <span style={{ flex: 1, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
          {job.label}
        </span>
        <span style={{ color: 'var(--text-muted)', flexShrink: 0, fontVariantNumeric: 'tabular-nums' }}>
          {duration}
        </span>
        {job.status === 'running' && (
          <button
            onClick={(e) => { e.stopPropagation(); cancel(job.id); }}
            style={{
              background: 'none', border: 'none', color: 'var(--accent-error, #ef4444)',
              cursor: 'pointer', fontSize: 11, padding: '0 2px', flexShrink: 0,
            }}
          >
            Cancel
          </button>
        )}
        {job.status === 'error' && job.error && (
          <span style={{ color: 'var(--accent-error, #ef4444)', flexShrink: 0 }} title={job.error}>
            err
          </span>
        )}
      </div>
      {/* Progress bar for running jobs with a known fraction */}
      {job.status === 'running' && job.progressFraction != null && (
        <div style={{ display: 'flex', alignItems: 'center', gap: 6, marginTop: 3, paddingLeft: 13 }}>
          <div style={{
            flex: 1, height: 3, borderRadius: 2,
            background: 'var(--bg-tertiary, #333)',
            overflow: 'hidden',
          }}>
            <div style={{
              height: '100%',
              width: `${Math.min(100, job.progressFraction * 100)}%`,
              background: 'var(--accent-primary, #3b82f6)',
              borderRadius: 2,
              transition: 'width 0.15s',
            }} />
          </div>
          {job.progress && (
            <span style={{ color: 'var(--text-muted)', fontSize: 10, flexShrink: 0, fontVariantNumeric: 'tabular-nums' }}>
              {job.progress}
            </span>
          )}
        </div>
      )}
    </div>
  );
}

export function SolverStatusIndicator() {
  const jobs = useSolverJobStore((s) => s.jobs);
  const clearHistory = useSolverJobStore((s) => s.clearHistory);
  const activeJob = useActiveJob();
  const [open, setOpen] = useState(false);
  const [, setNow] = useState(Date.now());
  const popoverRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!activeJob) return;
    const id = setInterval(() => setNow(Date.now()), 250);
    return () => clearInterval(id);
  }, [activeJob]);

  useEffect(() => {
    if (!open) return;
    const handler = (e: MouseEvent) => {
      if (popoverRef.current && !popoverRef.current.contains(e.target as Node)) {
        setOpen(false);
      }
    };
    window.addEventListener('mousedown', handler);
    return () => window.removeEventListener('mousedown', handler);
  }, [open]);

  const recentJobs = jobs.slice(-8).reverse();
  const runningJobs = jobs.filter(j => j.status === 'running');

  return (
    <div ref={popoverRef} style={{ position: 'relative' }}>
      {/* Popover */}
      {open && recentJobs.length > 0 && (
        <div style={{
          position: 'absolute', bottom: '100%', right: 0,
          width: 340,
          marginBottom: 6,
          background: 'var(--bg-primary)', border: '1px solid var(--border-color)',
          borderRadius: 6,
          padding: '8px 10px',
          maxHeight: 280, overflowY: 'auto',
          boxShadow: '0 -4px 16px rgba(0,0,0,0.15)',
          zIndex: 1000,
        }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 4 }}>
            <span style={{ fontSize: 11, fontWeight: 600, color: 'var(--text-secondary)' }}>
              Solver Queue {runningJobs.length > 0 && `(${runningJobs.length} running)`}
            </span>
            <button
              onClick={clearHistory}
              style={{
                background: 'none', border: 'none', color: 'var(--text-muted)',
                cursor: 'pointer', fontSize: 10,
              }}
            >
              Clear
            </button>
          </div>
          {recentJobs.map((job) => (
            <JobRow key={job.id} job={job} />
          ))}
        </div>
      )}

      {/* Inline indicator */}
      <button
        onClick={() => setOpen((o) => !o)}
        data-tour="solver-status"
        style={{
          display: 'inline-flex', alignItems: 'center', gap: 6,
          background: 'none', border: 'none',
          cursor: 'pointer',
          padding: '2px 6px',
          borderRadius: 4,
          fontSize: 11,
          color: 'var(--text-muted)',
          transition: 'background 0.15s',
        }}
        onMouseEnter={(e) => { e.currentTarget.style.background = 'rgba(0,100,60,0.08)'; }}
        onMouseLeave={(e) => { e.currentTarget.style.background = 'none'; }}
        title={activeJob ? `Running: ${activeJob.label}` : 'Solver ready'}
      >
        <span style={{
          width: 7, height: 7, borderRadius: '50%', flexShrink: 0,
          background: activeJob
            ? 'var(--accent-warning, #eab308)'
            : 'var(--accent-success, #22c55e)',
          animation: activeJob ? 'solver-pulse 1.2s ease-in-out infinite' : 'none',
        }} />
        {activeJob ? (
          <>
            <span style={{ maxWidth: 140, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
              {activeJob.label}
            </span>
            {activeJob.progressFraction != null && (
              <MiniProgressBar fraction={activeJob.progressFraction} />
            )}
            {activeJob.progress && (
              <span style={{ fontVariantNumeric: 'tabular-nums', fontSize: 10 }}>
                {activeJob.progress}
              </span>
            )}
          </>
        ) : (
          <span>Ready</span>
        )}
      </button>

      <style>{`
        @keyframes solver-pulse {
          0%, 100% { opacity: 1; }
          50% { opacity: 0.4; }
        }
      `}</style>
    </div>
  );
}
