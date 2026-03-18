/**
 * SolverStatusIndicator - Compact status dot+label that sits in the brand footer.
 * Click to expand a popover with recent job history and cancel buttons.
 */

import { useState, useEffect, useRef } from 'react';
import { useSolverJobStore, useActiveJob, type SolverJob } from '../stores/solverJobStore';

function elapsed(ms: number): string {
  if (ms < 1000) return `${ms}ms`;
  const s = ms / 1000;
  if (s < 60) return `${s.toFixed(1)}s`;
  return `${Math.floor(s / 60)}m ${Math.round(s % 60)}s`;
}

function JobRow({ job }: { job: SolverJob }) {
  const cancel = useSolverJobStore((s) => s.cancel);
  const duration = job.finishedAt
    ? elapsed(job.finishedAt - job.startedAt)
    : elapsed(Date.now() - job.startedAt);

  return (
    <div style={{
      display: 'flex', alignItems: 'center', gap: 8,
      padding: '4px 0',
      borderBottom: '1px solid var(--border-color)',
      fontSize: 11,
    }}>
      <span style={{
        width: 7, height: 7, borderRadius: '50%', flexShrink: 0,
        background:
          job.status === 'running' ? 'var(--accent-warning, #eab308)' :
          job.status === 'done' ? 'var(--accent-success, #22c55e)' :
          job.status === 'error' ? 'var(--accent-error, #ef4444)' :
          'var(--text-muted)',
      }} />
      <span style={{ flex: 1, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
        {job.label}
      </span>
      {job.progress && job.status === 'running' && (
        <span style={{ color: 'var(--text-muted)', flexShrink: 0 }}>{job.progress}</span>
      )}
      <span style={{ color: 'var(--text-muted)', flexShrink: 0 }}>{duration}</span>
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
  );
}

export function SolverStatusIndicator() {
  const jobs = useSolverJobStore((s) => s.jobs);
  const clearHistory = useSolverJobStore((s) => s.clearHistory);
  const activeJob = useActiveJob();
  const [open, setOpen] = useState(false);
  const [now, setNow] = useState(Date.now());
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

  const recentJobs = jobs.slice(-5).reverse();

  return (
    <div ref={popoverRef} style={{ position: 'relative' }}>
      {/* Popover */}
      {open && recentJobs.length > 0 && (
        <div style={{
          position: 'absolute', bottom: '100%', right: 0,
          width: 320,
          marginBottom: 6,
          background: 'var(--bg-primary)', border: '1px solid var(--border-color)',
          borderRadius: 6,
          padding: '8px 10px',
          maxHeight: 200, overflowY: 'auto',
          boxShadow: '0 -4px 16px rgba(0,0,0,0.15)',
          zIndex: 1000,
        }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 4 }}>
            <span style={{ fontSize: 11, fontWeight: 600, color: 'var(--text-secondary)' }}>Recent Jobs</span>
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
            <span style={{ maxWidth: 120, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
              {activeJob.label}
            </span>
            <span style={{ fontVariantNumeric: 'tabular-nums' }}>
              {elapsed(now - activeJob.startedAt)}
            </span>
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
