/**
 * MenuBar - Top menu bar with File/Edit/Window/Help menus
 */

import { useState, useRef, useEffect } from 'react';
import { DarkModeToggle } from './DarkModeToggle';
import { useUndoRedo } from '../hooks/useUndoRedo';
import { useOnboarding } from '../onboarding';

interface PanelInfo {
  id: string;
  name: string;
}

interface MenuBarProps {
  panels: PanelInfo[];
  closedPanels: Set<string>;
  onRestorePanel: (panelId: string) => void;
  onResetLayout: () => void;
  wasmStatus: 'loading' | 'ready' | 'error';
}

export function MenuBar({
  panels,
  closedPanels,
  onRestorePanel,
  onResetLayout,
  wasmStatus,
}: MenuBarProps) {
  const [activeMenu, setActiveMenu] = useState<string | null>(null);
  const menuRef = useRef<HTMLDivElement>(null);
  
  // Undo/redo functionality
  const { undo, redo, canUndo, canRedo } = useUndoRedo();
  
  // Onboarding
  const { startTour, resetAllTours } = useOnboarding();

  // Close menu when clicking outside
  useEffect(() => {
    const handleClickOutside = (e: MouseEvent) => {
      if (menuRef.current && !menuRef.current.contains(e.target as Node)) {
        setActiveMenu(null);
      }
    };
    document.addEventListener('mousedown', handleClickOutside);
    return () => document.removeEventListener('mousedown', handleClickOutside);
  }, []);

  const toggleMenu = (menuName: string) => {
    setActiveMenu(activeMenu === menuName ? null : menuName);
  };

  const isPanelVisible = (panelId: string) => !closedPanels.has(panelId);

  const handleTogglePanel = (panelId: string) => {
    if (closedPanels.has(panelId)) {
      onRestorePanel(panelId);
    }
    setActiveMenu(null);
  };

  const handleResetLayout = () => {
    onResetLayout();
    setActiveMenu(null);
  };

  return (
    <div
      ref={menuRef}
      style={{
        display: 'flex',
        alignItems: 'center',
        height: '36px',
        background: 'var(--bg-secondary)',
        borderBottom: '1px solid var(--border-color)',
        userSelect: 'none',
      }}
    >
      {/* App title */}
      <div
        style={{
          padding: '0 16px',
          fontWeight: 700,
          fontSize: '15px',
          background: 'linear-gradient(135deg, var(--accent-primary), var(--accent-secondary))',
          WebkitBackgroundClip: 'text',
          WebkitTextFillColor: 'transparent',
          backgroundClip: 'text',
          borderRight: '1px solid var(--border-color)',
          height: '100%',
          display: 'flex',
          alignItems: 'center',
        }}
      >
        FlexFoil
      </div>

      {/* File Menu */}
      <MenuButton
        label="File"
        isActive={activeMenu === 'file'}
        onClick={() => toggleMenu('file')}
      />
      {activeMenu === 'file' && (
        <MenuDropdown>
          <MenuItem label="New NACA..." disabled />
          <MenuItem label="Import .dat..." disabled />
          <MenuDivider />
          <MenuItem label="Export .dat..." disabled />
          <MenuItem label="Export SVG..." disabled />
        </MenuDropdown>
      )}

      {/* Edit Menu */}
      <MenuButton
        label="Edit"
        isActive={activeMenu === 'edit'}
        onClick={() => toggleMenu('edit')}
        dataTour="menu-edit"
      />
      {activeMenu === 'edit' && (
        <MenuDropdown style={{ left: '120px' }}>
          <MenuItem 
            label="Undo" 
            shortcut="Cmd+Z" 
            disabled={!canUndo}
            onClick={() => {
              undo();
              setActiveMenu(null);
            }}
          />
          <MenuItem 
            label="Redo" 
            shortcut="Cmd+Shift+Z" 
            disabled={!canRedo}
            onClick={() => {
              redo();
              setActiveMenu(null);
            }}
          />
        </MenuDropdown>
      )}

      {/* Window Menu */}
      <MenuButton
        label="Window"
        isActive={activeMenu === 'window'}
        onClick={() => toggleMenu('window')}
      />
      {activeMenu === 'window' && (
        <MenuDropdown style={{ left: '168px' }}>
          <div
            style={{
              padding: '4px 12px',
              fontSize: '11px',
              fontWeight: 600,
              color: 'var(--text-muted)',
              textTransform: 'uppercase',
              letterSpacing: '0.5px',
            }}
          >
            Panels
          </div>
          {panels.map((panel) => (
            <MenuItem
              key={panel.id}
              label={panel.name}
              checked={isPanelVisible(panel.id)}
              onClick={() => handleTogglePanel(panel.id)}
              disabled={isPanelVisible(panel.id)}
            />
          ))}
          <MenuDivider />
          <MenuItem label="Reset Layout" onClick={handleResetLayout} />
        </MenuDropdown>
      )}

      {/* Help Menu */}
      <MenuButton
        label="Help"
        isActive={activeMenu === 'help'}
        onClick={() => toggleMenu('help')}
        dataTour="menu-help"
      />
      {activeMenu === 'help' && (
        <MenuDropdown style={{ left: '232px' }}>
          <div
            style={{
              padding: '4px 12px',
              fontSize: '11px',
              fontWeight: 600,
              color: 'var(--text-muted)',
              textTransform: 'uppercase',
              letterSpacing: '0.5px',
            }}
          >
            Tutorials
          </div>
          <MenuItem
            label="Welcome Tour"
            onClick={() => {
              startTour('welcome');
              setActiveMenu(null);
            }}
          />
          <MenuItem
            label="Airfoil Editing Guide"
            onClick={() => {
              startTour('airfoilEditing');
              setActiveMenu(null);
            }}
          />
          <MenuItem
            label="Solving Guide"
            onClick={() => {
              startTour('solving');
              setActiveMenu(null);
            }}
          />
          <MenuDivider />
          <MenuItem
            label="Reset Tutorial Progress"
            onClick={() => {
              resetAllTours();
              setActiveMenu(null);
            }}
          />
          <MenuDivider />
          <MenuItem label="About FlexFoil" disabled />
        </MenuDropdown>
      )}

      {/* Spacer */}
      <div style={{ flex: 1 }} />

      {/* WASM Status */}
      <div
        style={{
          padding: '0 12px',
          fontSize: '11px',
          fontFamily: 'var(--font-mono)',
          color:
            wasmStatus === 'ready'
              ? 'var(--accent-primary)'
              : wasmStatus === 'error'
              ? 'var(--accent-danger)'
              : 'var(--accent-warning)',
          borderLeft: '1px solid var(--border-color)',
          height: '100%',
          display: 'flex',
          alignItems: 'center',
        }}
      >
        {wasmStatus === 'loading' && 'Loading...'}
        {wasmStatus === 'ready' && 'WASM Ready'}
        {wasmStatus === 'error' && 'WASM Error'}
      </div>

      {/* Dark Mode Toggle */}
      <div
        style={{
          padding: '0 8px',
          borderLeft: '1px solid var(--border-color)',
          height: '100%',
          display: 'flex',
          alignItems: 'center',
        }}
      >
        <DarkModeToggle />
      </div>
    </div>
  );
}

// Menu button component
function MenuButton({
  label,
  isActive,
  onClick,
  dataTour,
}: {
  label: string;
  isActive: boolean;
  onClick: () => void;
  dataTour?: string;
}) {
  const [isHovered, setIsHovered] = useState(false);

  return (
    <button
      onClick={onClick}
      onMouseEnter={() => setIsHovered(true)}
      onMouseLeave={() => setIsHovered(false)}
      data-tour={dataTour}
      style={{
        padding: '0 12px',
        height: '100%',
        background: isActive || isHovered ? 'var(--bg-hover)' : 'transparent',
        border: 'none',
        color: 'var(--text-secondary)',
        fontSize: '13px',
        cursor: 'pointer',
        transition: 'background 0.15s ease',
      }}
    >
      {label}
    </button>
  );
}

// Menu dropdown container
function MenuDropdown({
  children,
  style,
}: {
  children: React.ReactNode;
  style?: React.CSSProperties;
}) {
  return (
    <div
      style={{
        position: 'absolute',
        top: '36px',
        left: '80px',
        minWidth: '180px',
        background: 'var(--bg-secondary)',
        border: '1px solid var(--border-color)',
        borderRadius: '6px',
        boxShadow: '0 8px 24px rgba(0, 0, 0, 0.3)',
        zIndex: 1000,
        padding: '4px 0',
        ...style,
      }}
    >
      {children}
    </div>
  );
}

// Menu item component
function MenuItem({
  label,
  shortcut,
  checked,
  disabled,
  onClick,
}: {
  label: string;
  shortcut?: string;
  checked?: boolean;
  disabled?: boolean;
  onClick?: () => void;
}) {
  const [isHovered, setIsHovered] = useState(false);

  return (
    <button
      onClick={disabled ? undefined : onClick}
      onMouseEnter={() => setIsHovered(true)}
      onMouseLeave={() => setIsHovered(false)}
      disabled={disabled}
      style={{
        display: 'flex',
        alignItems: 'center',
        gap: '8px',
        width: '100%',
        padding: '6px 12px',
        background: isHovered && !disabled ? 'var(--bg-hover)' : 'transparent',
        border: 'none',
        color: disabled ? 'var(--text-muted)' : 'var(--text-secondary)',
        fontSize: '13px',
        textAlign: 'left',
        cursor: disabled ? 'default' : 'pointer',
        transition: 'background 0.15s ease',
      }}
    >
      {/* Checkmark */}
      <span style={{ width: '16px', color: 'var(--accent-primary)' }}>
        {checked !== undefined ? (checked ? '✓' : '') : ''}
      </span>

      {/* Label */}
      <span style={{ flex: 1 }}>{label}</span>

      {/* Shortcut */}
      {shortcut && (
        <span style={{ fontSize: '11px', color: 'var(--text-muted)' }}>{shortcut}</span>
      )}
    </button>
  );
}

// Menu divider
function MenuDivider() {
  return (
    <div
      style={{
        height: '1px',
        background: 'var(--border-color)',
        margin: '4px 0',
      }}
    />
  );
}
