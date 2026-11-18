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

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=assets/osabt.ico");

    if env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("windows") {
        return;
    }

    let manifest_dir = match env::var("CARGO_MANIFEST_DIR") {
        Ok(d) => PathBuf::from(d),
        Err(_) => return,
    };
    let out_dir = match env::var("OUT_DIR") {
        Ok(d) => PathBuf::from(d),
        Err(_) => return,
    };
    let ico_path = manifest_dir.join("assets").join("osabt.ico");
    if !ico_path.exists() {
        println!("cargo:warning=Missing assets/osabt.ico; skipping Windows PE resource embedding.");
        return;
    }

    let ico_rc = ico_path.to_string_lossy().replace('\\', "/");
    let rc_content = format!(
        "#pragma code_page(65001)\n\
         1 ICON \"{ico_rc}\"\n\
         1 VERSIONINFO\n\
         FILEVERSION 2,0,0,0\n\
         PRODUCTVERSION 2,0,0,0\n\
         FILEFLAGSMASK 0x3fL\n\
         FILEFLAGS 0x0L\n\
         FILEOS 0x40004L\n\
         FILETYPE 0x1L\n\
         FILESUBTYPE 0x0L\n\
         BEGIN\n\
           BLOCK \"StringFileInfo\"\n\
           BEGIN\n\
             BLOCK \"040904b0\"\n\
             BEGIN\n\
               VALUE \"CompanyName\", \"CRV Science Network\\0\"\n\
               VALUE \"FileDescription\", \"OSABT — Open Source Ad-hoc Biology Tools\\0\"\n\
               VALUE \"FileVersion\", \"2.0.0.0\\0\"\n\
               VALUE \"InternalName\", \"osabt\\0\"\n\
               VALUE \"LegalCopyright\", \"Copyright CRV Science Network (GPL-3.0-or-later)\\0\"\n\
               VALUE \"OriginalFilename\", \"osabt.exe\\0\"\n\
               VALUE \"ProductName\", \"OSABT\\0\"\n\
               VALUE \"ProductVersion\", \"2.0.0\\0\"\n\
             END\n\
           END\n\
           BLOCK \"VarFileInfo\"\n\
           BEGIN\n\
             VALUE \"Translation\", 0x0409, 1200\n\
           END\n\
         END\n"
    );

    let rc_path = out_dir.join("osabt_res.rc");
    if let Err(err) = fs::write(&rc_path, rc_content) {
        println!("cargo:warning=Failed writing {}: {err}", rc_path.display());
        return;
    }

    let target_env = env::var("CARGO_CFG_TARGET_ENV").unwrap_or_default();
    if target_env == "gnu" {
        let obj_path = out_dir.join("osabt_res.o");
        let windres_bin = ["windres", "x86_64-w64-mingw32-windres"]
            .into_iter()
            .find(|bin| Command::new(bin).arg("--version").output().map(|o| o.status.success()).unwrap_or(false));
        if let Some(bin) = windres_bin {
            match Command::new(bin)
                .arg("-i")
                .arg(&rc_path)
                .arg("-o")
                .arg(&obj_path)
                .arg("-O")
                .arg("coff")
                .output()
            {
                Ok(out) if out.status.success() && obj_path.exists() => {
                    println!("cargo:rustc-link-arg-bins={}", obj_path.display());
                }
                Ok(out) => {
                    println!(
                        "cargo:warning={bin} failed: {}",
                        String::from_utf8_lossy(&out.stderr).trim()
                    );
                }
                Err(err) => println!("cargo:warning=Failed running {bin}: {err}"),
            }
        } else {
            println!("cargo:warning=windres not found; skipping Windows PE icon/version resources.");
        }
    } else if target_env == "msvc" {
        let res_path = out_dir.join("osabt_res.res");
        if let Some(rc_bin) = find_msvc_rc() {
            match Command::new(&rc_bin)
                .arg("/nologo")
                .arg("/fo")
                .arg(&res_path)
                .arg(&rc_path)
                .output()
            {
                Ok(out) if out.status.success() && res_path.exists() => {
                    println!("cargo:rustc-link-arg-bins={}", res_path.display());
                }
                Ok(out) => {
                    println!(
                        "cargo:warning=rc.exe failed: {}",
                        String::from_utf8_lossy(&out.stderr).trim()
                    );
                }
                Err(err) => println!("cargo:warning=Failed running rc.exe: {err}"),
            }
        } else {
            println!("cargo:warning=rc.exe not found; skipping Windows PE icon/version resources.");
        }
    }
}

fn find_msvc_rc() -> Option<PathBuf> {
    if Command::new("rc.exe")
        .arg("/?")
        .output()
        .map(|o| o.status.success() || !o.stdout.is_empty() || !o.stderr.is_empty())
        .unwrap_or(false)
    {
        return Some(PathBuf::from("rc.exe"));
    }
    if let Ok(sdk_dir) = env::var("WindowsSdkDir") {
        let bin_root = Path::new(&sdk_dir).join("bin");
        if let Some(found) = scan_windows_kits_rc(&bin_root) {
            return Some(found);
        }
    }
    for base in [
        r"C:\Program Files (x86)\Windows Kits\10\bin",
        r"C:\Program Files\Windows Kits\10\bin",
    ] {
        if let Some(found) = scan_windows_kits_rc(Path::new(base)) {
            return Some(found);
        }
    }
    None
}

fn scan_windows_kits_rc(bin_root: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(bin_root).ok()?;
    let mut candidates: Vec<PathBuf> = entries
        .filter_map(|e| e.ok().map(|x| x.path()))
        .filter(|p| p.is_dir())
        .collect();
    candidates.sort();
    candidates.reverse();
    for dir in candidates {
        for arch in ["x64", "x86", "arm64"] {
            let candidate = dir.join(arch).join("rc.exe");
            if candidate.exists() {
                return Some(candidate);
            }
        }
    }
    None
}
