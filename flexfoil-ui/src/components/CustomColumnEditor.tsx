import { useState, useCallback, useMemo } from 'react';
import {
  useCustomColumnStore,
  type CustomColumn,
} from '../stores/customColumnStore';
import {
  validateExpression,
  EXPRESSION_COLUMNS,
  AVAILABLE_FUNCTIONS,
} from '../lib/expressionEngine';

interface CustomColumnEditorProps {
  onClose: () => void;
}

export function CustomColumnEditor({ onClose }: CustomColumnEditorProps) {
  const { columns, addColumn, removeColumn } = useCustomColumnStore();
  const [name, setName] = useState('');
  const [expression, setExpression] = useState('');
  const [error, setError] = useState<string | null>(null);
  const [confirmDeleteId, setConfirmDeleteId] = useState<string | null>(null);

  const validation = useMemo(() => {
    if (!expression.trim()) return null;
    return validateExpression(expression);
  }, [expression]);

  const handleAdd = useCallback(() => {
    const trimmedName = name.trim();
    if (!trimmedName) {
      setError('Name is required');
      return;
    }
    if (!expression.trim()) {
      setError('Expression is required');
      return;
    }
    const err = addColumn(trimmedName, expression);
    if (err) {
      setError(err);
      return;
    }
    setName('');
    setExpression('');
    setError(null);
  }, [name, expression, addColumn]);

  const handleDelete = useCallback(
    (id: string) => {
      if (confirmDeleteId === id) {
        removeColumn(id);
        setConfirmDeleteId(null);
      } else {
        setConfirmDeleteId(id);
      }
    },
    [confirmDeleteId, removeColumn],
  );

  const inputStyle: React.CSSProperties = {
    padding: '5px 8px',
    fontSize: '12px',
    background: 'var(--bg-tertiary)',
    border: '1px solid var(--border-color)',
    borderRadius: '4px',
    color: 'var(--text-primary)',
    width: '100%',
    boxSizing: 'border-box',
  };

  return (
    <div
      onClick={onClose}
      style={{
        position: 'fixed',
        inset: 0,
        background: 'rgba(0,0,0,0.4)',
        zIndex: 9999,
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
      }}
    >
      <div
        onClick={(e) => e.stopPropagation()}
        style={{
          background: 'var(--bg-secondary)',
          border: '1px solid var(--border-color)',
          borderRadius: '10px',
          boxShadow: '0 20px 60px rgba(0,0,0,0.3)',
          width: '480px',
          maxWidth: '90vw',
          maxHeight: '80vh',
          display: 'flex',
          flexDirection: 'column',
        }}
      >
        {/* Header */}
        <div
          style={{
            padding: '16px 20px 12px',
            borderBottom: '1px solid var(--border-color)',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            flexShrink: 0,
          }}
        >
          <div
            style={{
              fontSize: '14px',
              fontWeight: 700,
              color: 'var(--text-primary)',
            }}
          >
            Custom Columns
          </div>
          <button
            onClick={onClose}
            style={{
              background: 'none',
              border: 'none',
              fontSize: '16px',
              color: 'var(--text-muted)',
              cursor: 'pointer',
              padding: '2px 6px',
            }}
          >
            &times;
          </button>
        </div>

        {/* Content */}
        <div
          style={{
            flex: 1,
            overflow: 'auto',
            padding: '16px 20px',
            display: 'flex',
            flexDirection: 'column',
            gap: '14px',
          }}
        >
          {/* New column form */}
          <div
            style={{
              display: 'flex',
              flexDirection: 'column',
              gap: '8px',
              padding: '12px',
              background: 'var(--bg-primary)',
              borderRadius: '8px',
              border: '1px solid var(--border-color)',
            }}
          >
            <div style={{ fontSize: '11px', fontWeight: 600, color: 'var(--text-muted)', textTransform: 'uppercase', letterSpacing: '0.5px' }}>
              Add Column
            </div>

            <div style={{ display: 'flex', gap: '8px' }}>
              <div style={{ flex: '0 0 130px' }}>
                <input
                  placeholder="Column name"
                  value={name}
                  onChange={(e) => { setName(e.target.value); setError(null); }}
                  onKeyDown={(e) => { if (e.key === 'Enter') handleAdd(); }}
                  style={inputStyle}
                />
              </div>
              <div style={{ flex: 1 }}>
                <input
                  placeholder="e.g. cl^2 / cd"
                  value={expression}
                  onChange={(e) => { setExpression(e.target.value); setError(null); }}
                  onKeyDown={(e) => { if (e.key === 'Enter') handleAdd(); }}
                  style={{
                    ...inputStyle,
                    fontFamily: 'var(--font-mono)',
                    borderColor:
                      expression.trim()
                        ? validation?.valid
                          ? 'var(--accent-primary)'
                          : 'var(--accent-danger)'
                        : 'var(--border-color)',
                  }}
                />
              </div>
            </div>

            {/* Validation / error feedback */}
            {(error || (expression.trim() && validation && !validation.valid)) && (
              <div style={{ fontSize: '11px', color: 'var(--accent-danger)' }}>
                {error || validation?.error}
              </div>
            )}

            <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
              <button
                onClick={handleAdd}
                disabled={!name.trim() || !expression.trim() || (validation != null && !validation.valid)}
                style={{
                  padding: '5px 14px',
                  fontSize: '11px',
                  fontWeight: 600,
                  background: 'var(--accent-primary)',
                  border: 'none',
                  borderRadius: '5px',
                  color: 'var(--bg-primary)',
                  cursor: 'pointer',
                  opacity: (!name.trim() || !expression.trim() || (validation != null && !validation.valid)) ? 0.4 : 1,
                }}
              >
                Add
              </button>
              {expression.trim() && validation?.valid && (
                <span style={{ fontSize: '10px', color: 'var(--accent-primary)' }}>
                  Valid expression
                </span>
              )}
            </div>

            {/* Reference */}
            <details style={{ fontSize: '10px', color: 'var(--text-muted)' }}>
              <summary style={{ cursor: 'pointer', userSelect: 'none' }}>
                Available columns &amp; functions
              </summary>
              <div style={{ marginTop: '6px', lineHeight: 1.7 }}>
                <div>
                  <strong>Columns:</strong>{' '}
                  <span style={{ fontFamily: 'var(--font-mono)' }}>
                    {EXPRESSION_COLUMNS.join(', ')}
                  </span>
                </div>
                <div>
                  <strong>Functions:</strong>{' '}
                  <span style={{ fontFamily: 'var(--font-mono)' }}>
                    {AVAILABLE_FUNCTIONS.join(', ')}
                  </span>
                </div>
                <div>
                  <strong>Operators:</strong>{' '}
                  <span style={{ fontFamily: 'var(--font-mono)' }}>
                    + - * / ^ % ( )
                  </span>
                </div>
                <div>
                  <strong>Constants:</strong>{' '}
                  <span style={{ fontFamily: 'var(--font-mono)' }}>pi, e</span>
                </div>
              </div>
            </details>
          </div>

          {/* Existing columns */}
          {columns.length > 0 && (
            <div style={{ display: 'flex', flexDirection: 'column', gap: '6px' }}>
              <div
                style={{
                  fontSize: '11px',
                  fontWeight: 600,
                  color: 'var(--text-muted)',
                  textTransform: 'uppercase',
                  letterSpacing: '0.5px',
                }}
              >
                Active Columns ({columns.length})
              </div>
              {columns.map((col) => (
                <CustomColumnRow
                  key={col.id}
                  column={col}
                  isConfirmingDelete={confirmDeleteId === col.id}
                  onDelete={() => handleDelete(col.id)}
                  onCancelDelete={() => setConfirmDeleteId(null)}
                />
              ))}
            </div>
          )}

          {columns.length === 0 && (
            <div
              style={{
                textAlign: 'center',
                padding: '20px',
                color: 'var(--text-muted)',
                fontSize: '12px',
              }}
            >
              No custom columns yet. Add one above.
            </div>
          )}
        </div>
      </div>
    </div>
  );
}

function CustomColumnRow({
  column,
  isConfirmingDelete,
  onDelete,
  onCancelDelete,
}: {
  column: CustomColumn;
  isConfirmingDelete: boolean;
  onDelete: () => void;
  onCancelDelete: () => void;
}) {
  return (
    <div
      style={{
        display: 'flex',
        alignItems: 'center',
        gap: '8px',
        padding: '8px 10px',
        background: 'var(--bg-primary)',
        borderRadius: '6px',
        border: '1px solid var(--border-color)',
      }}
    >
      <div style={{ flex: 1, minWidth: 0 }}>
        <div
          style={{
            fontSize: '12px',
            fontWeight: 600,
            color: 'var(--text-primary)',
            overflow: 'hidden',
            textOverflow: 'ellipsis',
            whiteSpace: 'nowrap',
          }}
        >
          {column.name}
        </div>
        <div
          style={{
            fontSize: '11px',
            fontFamily: 'var(--font-mono)',
            color: 'var(--text-muted)',
            overflow: 'hidden',
            textOverflow: 'ellipsis',
            whiteSpace: 'nowrap',
          }}
        >
          {column.expression}
        </div>
      </div>

      {isConfirmingDelete ? (
        <div style={{ display: 'flex', gap: '4px', flexShrink: 0 }}>
          <button
            onClick={onDelete}
            style={{
              background: 'var(--accent-danger)',
              border: 'none',
              borderRadius: '3px',
              padding: '3px 8px',
              fontSize: '10px',
              color: '#fff',
              cursor: 'pointer',
            }}
          >
            Delete
          </button>
          <button
            onClick={onCancelDelete}
            style={{
              background: 'transparent',
              border: '1px solid var(--border-color)',
              borderRadius: '3px',
              padding: '3px 6px',
              fontSize: '10px',
              color: 'var(--text-muted)',
              cursor: 'pointer',
            }}
          >
            No
          </button>
        </div>
      ) : (
        <button
          onClick={onDelete}
          title="Remove column"
          style={{
            background: 'none',
            border: 'none',
            cursor: 'pointer',
            padding: '2px 4px',
            fontSize: '12px',
            opacity: 0.3,
            color: 'var(--text-primary)',
            flexShrink: 0,
          }}
          onMouseEnter={(e) => {
            e.currentTarget.style.opacity = '1';
          }}
          onMouseLeave={(e) => {
            e.currentTarget.style.opacity = '0.3';
          }}
        >
          &times;
        </button>
      )}
    </div>
  );
}
