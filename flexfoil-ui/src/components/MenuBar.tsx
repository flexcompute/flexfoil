/**
 * MenuBar - Top menu bar with File/Edit/Window/Help menus
 */

import { useState, useRef, useEffect } from 'react';
import { DarkModeToggle } from './DarkModeToggle';
import { useUndoRedo } from '../hooks/useUndoRedo';
import { useOnboarding } from '../onboarding';
import { FlexcomputeLogo } from './FlexcomputeLogo';
import { AboutDialog } from './AboutDialog';
import { ChangelogDialog } from './ChangelogDialog';

const DOCUMENTATION_URL = 'https://foil.flexcompute.com/docs/';
const FLEXCOMPUTE_URL = 'https://www.flexcompute.com/';

interface PanelInfo {
  id: string;
  name: string;
}

interface MenuBarProps {
  panels: PanelInfo[];
  closedPanels: Set<string>;
  onRestorePanel: (panelId: string) => void;
  onOpenPanel: (panelId: string) => void;  // Focus/select an existing panel
  onResetLayout: () => void;
  wasmStatus: 'loading' | 'ready' | 'error';
}

export function MenuBar({
  panels,
  closedPanels,
  onRestorePanel,
  onOpenPanel,
  onResetLayout,
  wasmStatus,
}: MenuBarProps) {
  const [activeMenu, setActiveMenu] = useState<string | null>(null);
  const [showAbout, setShowAbout] = useState(false);
  const [showChangelog, setShowChangelog] = useState(false);
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
      // Panel is closed - restore it
      onRestorePanel(panelId);
    } else {
      // Panel exists but might be hidden behind another tab - bring to front
      onOpenPanel(panelId);
    }
    setActiveMenu(null);
  };

  const handleResetLayout = () => {
    onResetLayout();
    setActiveMenu(null);
  };

  const openDocumentation = () => {
    window.open(DOCUMENTATION_URL, '_blank', 'noopener,noreferrer');
    setActiveMenu(null);
  };

  return (
    <div
      ref={menuRef}
      style={{
        display: 'flex',
        alignItems: 'center',
        height: '44px',
        background: 'var(--bg-secondary)',
        borderBottom: '1px solid var(--border-color)',
        userSelect: 'none',
      }}
    >
      {/* Brand */}
      <a
        href={FLEXCOMPUTE_URL}
        target="_blank"
        rel="noopener noreferrer"
        aria-label="Open Flexcompute.com"
        style={{
          display: 'flex',
          alignItems: 'center',
          gap: '6px',
          padding: '0 16px',
          borderRight: '1px solid var(--border-color)',
          height: '100%',
          minWidth: '220px',
          textDecoration: 'none',
          color: 'var(--text-primary)',
          background: 'var(--brand-header-surface)',
        }}
      >
        <span
          style={{
            fontSize: '13px',
            fontWeight: 700,
            letterSpacing: '0.08em',
            textTransform: 'uppercase',
            color: 'var(--brand-primary)',
          }}
        >
          FlexFoil
        </span>
        <span
          style={{
            fontSize: '10px',
            color: 'var(--text-secondary)',
            alignSelf: 'center',
          }}
        >
          by
        </span>
        <FlexcomputeLogo height={18} />
      </a>

      {/* File Menu */}
      <MenuGroup>
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
      </MenuGroup>

      {/* Edit Menu */}
      <MenuGroup>
        <MenuButton
          label="Edit"
          isActive={activeMenu === 'edit'}
          onClick={() => toggleMenu('edit')}
          dataTour="menu-edit"
        />
        {activeMenu === 'edit' && (
          <MenuDropdown>
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
      </MenuGroup>

      {/* Window Menu */}
      <MenuGroup>
        <MenuButton
          label="Window"
          isActive={activeMenu === 'window'}
          onClick={() => toggleMenu('window')}
          dataTour="menu-window"
        />
        {activeMenu === 'window' && (
          <MenuDropdown>
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
              />
            ))}
            <MenuDivider />
            <MenuItem label="Reset Layout" onClick={handleResetLayout} />
          </MenuDropdown>
        )}
      </MenuGroup>

      {/* Help Menu */}
      <MenuGroup>
        <MenuButton
          label="Help"
          isActive={activeMenu === 'help'}
          onClick={() => toggleMenu('help')}
          dataTour="menu-help"
        />
        {activeMenu === 'help' && (
          <MenuDropdown>
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
            <MenuItem label="Documentation" onClick={openDocumentation} />
            <MenuDivider />
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
            <MenuItem
              label="Data Explorer Guide"
              onClick={() => {
                startTour('dataExplorer', true);
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
            <MenuItem
              label="What's New"
              onClick={() => {
                setShowChangelog(true);
                setActiveMenu(null);
              }}
            />
            <MenuItem
              label="About FlexFoil"
              onClick={() => {
                setShowAbout(true);
                setActiveMenu(null);
              }}
            />
          </MenuDropdown>
        )}
      </MenuGroup>

      <MenuButton
        label="Docs"
        isActive={false}
        onClick={openDocumentation}
      />

      {/* Spacer */}
      <div style={{ flex: 1 }} />

      {/* Restore Tutorial Button - shows after tutorial has been dismissed */}
      <RestoreTutorialButton />

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

      <AboutDialog
        open={showAbout}
        onClose={() => setShowAbout(false)}
        onOpenChangelog={() => setShowChangelog(true)}
      />
      <ChangelogDialog
        open={showChangelog}
        onClose={() => setShowChangelog(false)}
      />
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
        fontWeight: 500,
        cursor: 'pointer',
        transition: 'background 0.15s ease',
      }}
    >
      {label}
    </button>
  );
}

function MenuGroup({ children }: { children: React.ReactNode }) {
  return (
    <div
      style={{
        position: 'relative',
        height: '100%',
        display: 'flex',
        alignItems: 'stretch',
      }}
    >
      {children}
    </div>
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
        top: '100%',
        left: 0,
        minWidth: '180px',
        background: 'var(--bg-secondary)',
        border: '1px solid var(--border-color)',
        borderRadius: '6px',
        boxShadow: '0 16px 40px rgba(0, 0, 0, 0.22)',
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

// Restore Tutorial button - appears after tutorial is dismissed
function RestoreTutorialButton() {
  const [isHovered, setIsHovered] = useState(false);
  const { hasStartedTour, startTour, isActive } = useOnboarding();
  
  // Only show if welcome tour has been started (dismissed or completed) and no tour is active
  if (!hasStartedTour('welcome') || isActive) {
    return null;
  }

  return (
    <button
      onClick={() => startTour('welcome')}
      onMouseEnter={() => setIsHovered(true)}
      onMouseLeave={() => setIsHovered(false)}
      title="Restart the welcome tutorial"
      style={{
        padding: '0 10px',
        height: '24px',
        marginRight: '8px',
        background: isHovered 
          ? 'var(--accent-primary)' 
          : 'transparent',
        border: '1px solid var(--accent-primary)',
        borderRadius: '4px',
        color: isHovered ? 'white' : 'var(--accent-primary)',
        fontSize: '11px',
        fontWeight: 500,
        cursor: 'pointer',
        transition: 'all 0.15s ease',
        display: 'flex',
        alignItems: 'center',
        gap: '4px',
      }}
    >
      <span style={{ fontSize: '12px' }}>?</span>
      Tutorial
    </button>
  );
}
