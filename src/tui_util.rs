// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Shared RAII Ratatui + Crossterm alternate-screen terminal lifecycle guard (`TerminalSession`)
//! and session runners (`run_tui_session`, `run_tui_session_io`) across OSABT TUI workbenches.

use std::io;

use crossterm::execute;
use crossterm::terminal::{
    disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen,
};
use ratatui::backend::CrosstermBackend;
use ratatui::Terminal;

/// Standard stdout-backed Ratatui terminal type used by OSABT interactive TUIs.
pub type StdoutTerminal = Terminal<CrosstermBackend<io::Stdout>>;

/// RAII terminal guard that enters raw mode and the alternate screen on construction
/// and restores terminal state on explicit completion or `Drop`.
pub struct TerminalSession {
    terminal: StdoutTerminal,
    cleaned_up: bool,
}

impl TerminalSession {
    /// Enters terminal raw mode and the alternate screen, returning `io::Result`.
    pub fn enter_io() -> io::Result<Self> {
        enable_raw_mode()?;
        let mut stdout = io::stdout();
        if let Err(err) = execute!(stdout, EnterAlternateScreen) {
            let _ = disable_raw_mode();
            return Err(err);
        }
        let backend = CrosstermBackend::new(stdout);
        match Terminal::new(backend) {
            Ok(terminal) => Ok(Self {
                terminal,
                cleaned_up: false,
            }),
            Err(err) => {
                let mut err_stdout = io::stdout();
                let _ = disable_raw_mode();
                let _ = execute!(err_stdout, LeaveAlternateScreen);
                Err(err)
            }
        }
    }

    /// Enters terminal raw mode and the alternate screen with descriptive `String` errors.
    pub fn enter() -> Result<Self, String> {
        enable_raw_mode().map_err(|e| format!("Failed enabling terminal raw mode: {e}"))?;
        let mut stdout = io::stdout();
        if let Err(e) = execute!(stdout, EnterAlternateScreen) {
            let _ = disable_raw_mode();
            return Err(format!("Failed entering alternate screen: {e}"));
        }
        let backend = CrosstermBackend::new(stdout);
        match Terminal::new(backend) {
            Ok(terminal) => Ok(Self {
                terminal,
                cleaned_up: false,
            }),
            Err(e) => {
                let mut err_stdout = io::stdout();
                let _ = disable_raw_mode();
                let _ = execute!(err_stdout, LeaveAlternateScreen);
                Err(format!("Failed creating Ratatui terminal: {e}"))
            }
        }
    }

    /// Returns a mutable borrow of the wrapped `StdoutTerminal`.
    pub fn terminal_mut(&mut self) -> &mut StdoutTerminal {
        &mut self.terminal
    }

    /// Leaves the alternate screen, disables raw mode, and restores the cursor,
    /// propagating any `io::Error` encountered after attempting all cleanup steps.
    pub fn finish_io(mut self) -> io::Result<()> {
        self.cleaned_up = true;
        let r1 = disable_raw_mode();
        let r2 = execute!(self.terminal.backend_mut(), LeaveAlternateScreen);
        let r3 = self.terminal.show_cursor();
        r1?;
        r2?;
        r3?;
        Ok(())
    }

    /// Best-effort quiet restoration of terminal state without propagating errors.
    pub fn finish_quiet(mut self) {
        self.restore_quiet();
    }

    fn restore_quiet(&mut self) {
        if !self.cleaned_up {
            self.cleaned_up = true;
            let _ = disable_raw_mode();
            let _ = execute!(self.terminal.backend_mut(), LeaveAlternateScreen);
            let _ = self.terminal.show_cursor();
        }
    }
}

impl Drop for TerminalSession {
    fn drop(&mut self) {
        self.restore_quiet();
    }
}

/// Runs an interactive TUI event loop inside a `TerminalSession` with `Result<T, String>` errors.
pub fn run_tui_session<T, F>(loop_fn: F) -> Result<T, String>
where
    F: FnOnce(&mut StdoutTerminal) -> Result<T, String>,
{
    let mut session = TerminalSession::enter()?;
    let res = loop_fn(session.terminal_mut());
    session.finish_quiet();
    res
}

/// Runs an interactive TUI event loop inside a `TerminalSession` with `io::Result<T>` errors.
pub fn run_tui_session_io<T, F>(loop_fn: F) -> io::Result<T>
where
    F: FnOnce(&mut StdoutTerminal) -> io::Result<T>,
{
    let mut session = TerminalSession::enter_io()?;
    let res = loop_fn(session.terminal_mut());
    session.finish_io()?;
    res
}
