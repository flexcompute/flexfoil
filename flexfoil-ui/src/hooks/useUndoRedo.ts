/**
 * useUndoRedo - Hook for undo/redo functionality with keyboard shortcuts
 */

import { useEffect, useCallback } from 'react';
import { 
  useTemporalStore, 
  undo as performUndo, 
  redo as performRedo 
} from '../stores/airfoilStore';

/**
 * Hook that provides undo/redo state and keyboard shortcuts.
 * 
 * Keyboard shortcuts:
 * - Cmd/Ctrl + Z: Undo
 * - Cmd/Ctrl + Shift + Z: Redo
 * - Cmd/Ctrl + Y: Redo (Windows style)
 */
export function useUndoRedo() {
  // Get temporal state
  const canUndo = useTemporalStore((state) => state.pastStates.length > 0);
  const canRedo = useTemporalStore((state) => state.futureStates.length > 0);
  const pastStatesCount = useTemporalStore((state) => state.pastStates.length);
  const futureStatesCount = useTemporalStore((state) => state.futureStates.length);

  // Wrapped undo/redo functions with checks
  const undo = useCallback(() => {
    if (canUndo) {
      performUndo();
    }
  }, [canUndo]);

  const redo = useCallback(() => {
    if (canRedo) {
      performRedo();
    }
  }, [canRedo]);

  // Keyboard shortcut handler
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      // Check if user is typing in an input/textarea
      const target = e.target as HTMLElement;
      if (
        target.tagName === 'INPUT' ||
        target.tagName === 'TEXTAREA' ||
        target.isContentEditable
      ) {
        return;
      }

      const isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
      const modifier = isMac ? e.metaKey : e.ctrlKey;

      if (!modifier) return;

      // Cmd/Ctrl + Z
      if (e.key === 'z' || e.key === 'Z') {
        e.preventDefault();
        if (e.shiftKey) {
          // Cmd/Ctrl + Shift + Z: Redo
          redo();
        } else {
          // Cmd/Ctrl + Z: Undo
          undo();
        }
        return;
      }

      // Cmd/Ctrl + Y: Redo (Windows style)
      if (e.key === 'y' || e.key === 'Y') {
        e.preventDefault();
        redo();
        return;
      }
    };

    window.addEventListener('keydown', handleKeyDown);
    return () => window.removeEventListener('keydown', handleKeyDown);
  }, [undo, redo]);

  return {
    undo,
    redo,
    canUndo,
    canRedo,
    pastStatesCount,
    futureStatesCount,
  };
}

export default useUndoRedo;
