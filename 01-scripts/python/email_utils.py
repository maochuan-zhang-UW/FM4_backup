"""
email_utils.py
==============
Email notifications for the FM7 pipeline via macOS Mail.app (osascript).
No SMTP credentials needed — uses whatever account is configured in Mail.app.
Pattern adapted from /Users/mczhang/Documents/GitHub/Axial_Redpy/scripts/setup/send_progress_email.py
"""
import subprocess
import tempfile
import os
from datetime import datetime
from pathlib import Path


def send_email(to: str, subject: str, body: str,
               attachments: list[str | Path] | None = None) -> bool:
    """
    Send an email via macOS Mail.app using AppleScript.

    Parameters
    ----------
    to          : recipient email address
    subject     : email subject line
    body        : plain-text body
    attachments : list of file paths (PNG plots etc.) to attach

    Returns
    -------
    True on success, False on failure (also prints the error).
    """
    # Escape characters that break AppleScript string literals
    def _esc(s: str) -> str:
        return s.replace('\\', '\\\\').replace('"', '\\"')

    body_esc    = _esc(body)
    subject_esc = _esc(subject)

    # Build attachment lines
    attach_lines = ''
    if attachments:
        for path in attachments:
            path = str(Path(path).resolve())
            attach_lines += f'\n        make new attachment with properties {{file name:(POSIX file "{path}" as alias)}}'

    script = f'''tell application "Mail"
    set theMessage to make new outgoing message with properties {{subject:"{subject_esc}", content:"{body_esc}", visible:false}}
    tell theMessage
        make new to recipient with properties {{address:"{to}"}}{attach_lines}
    end tell
    send theMessage
end tell'''

    with tempfile.NamedTemporaryFile(suffix='.applescript', mode='w', delete=False) as f:
        f.write(script)
        tmp = f.name
    try:
        subprocess.run(['osascript', tmp], check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f'[email_utils] Failed to send email: {e.stderr.decode()}')
        return False
    finally:
        os.unlink(tmp)


def step_email(step_num: int, step_name: str,
               stats: dict, plots: list[str | Path] | None = None,
               review_required: bool = False,
               error: str | None = None):
    """
    Send a step-completion email with stats and optional plot attachments.

    Parameters
    ----------
    step_num         : notebook number (1–12)
    step_name        : human-readable name
    stats            : dict of key→value pairs shown in the email body
    plots            : list of PNG file paths to attach
    review_required  : if True, adds a prominent review prompt in the subject
    error            : if not None, marks the email as an error report
    """
    from pipeline_config import EMAIL_TO

    now = datetime.now().strftime('%Y-%m-%d %H:%M')

    if error:
        subject = f'[FM7 ❌ ERROR] Step {step_num:02d} — {step_name}'
    elif review_required:
        subject = f'[FM7 ⚠️ REVIEW NEEDED] Step {step_num:02d} — {step_name}'
    else:
        subject = f'[FM7 ✓ Done] Step {step_num:02d}/{12} — {step_name}'

    lines = [
        f'FM7 Y1+Y2 Pipeline — Step {step_num:02d}: {step_name}',
        '=' * 55,
        f'Completed at: {now}',
        '',
    ]

    if error:
        lines += ['ERROR:', error, '']

    if stats:
        lines.append('--- Results ---')
        for k, v in stats.items():
            lines.append(f'  {k:<35} {v}')
        lines.append('')

    if review_required:
        lines += [
            '*** ACTION REQUIRED ***',
            'Please review the attached plots.',
            'When satisfied, set "awaiting_review": false in:',
            '  02-data/pipeline_progress.json',
            'Then re-run the pipeline to continue.',
            '',
        ]

    if plots:
        lines.append(f'Attached: {len(plots)} plot(s)')

    body = '\n'.join(lines)
    send_email(EMAIL_TO, subject, body, attachments=plots)
    print(f'[email_utils] Step email sent: {subject}')


def hourly_email(current_step: int, current_step_name: str,
                 completed_steps: list[dict],
                 progress_pct: float | None,
                 last_log_lines: list[str],
                 plots: list[str | Path] | None = None):
    """
    Send the hourly status update email.

    Parameters
    ----------
    current_step       : notebook number currently running
    current_step_name  : name of current notebook
    completed_steps    : list of {'step', 'name', 'duration_s', 'summary'} dicts
    progress_pct       : 0–100, current notebook progress (None if unknown)
    last_log_lines     : recent log lines to include
    plots              : optional plots from current step
    """
    from pipeline_config import EMAIL_TO

    now = datetime.now().strftime('%Y-%m-%d %H:%M')
    pct_str = f'{progress_pct:.0f}%' if progress_pct is not None else '?'

    subject = (f'[FM7 ↻ Hourly] Step {current_step:02d} — '
               f'{current_step_name} ({pct_str})')

    lines = [
        'FM7 Y1+Y2 Pipeline — Hourly Update',
        '=' * 55,
        f'Time: {now}',
        f'Current step: {current_step:02d} — {current_step_name}  ({pct_str} done)',
        '',
        '--- Completed steps ---',
    ]
    for s in completed_steps:
        dur = f"{s.get('duration_s', 0)/60:.1f} min"
        lines.append(f"  ✓ {s['step']:02d} {s['name']:<35} ({dur})  {s.get('summary','')}")

    if last_log_lines:
        lines += ['', '--- Recent log ---']
        lines += [f'  {l}' for l in last_log_lines[-15:]]

    body = '\n'.join(lines)
    send_email(EMAIL_TO, subject, body, attachments=plots)
    print(f'[email_utils] Hourly email sent: {subject}')


def error_email(step_num: int, step_name: str, traceback_str: str):
    """Send an immediate error alert — pipeline has stopped."""
    from pipeline_config import EMAIL_TO

    subject = f'[FM7 ❌ STOPPED] Step {step_num:02d} — {step_name} FAILED'
    body = '\n'.join([
        f'FM7 pipeline stopped at step {step_num:02d}: {step_name}',
        f'Time: {datetime.now().strftime("%Y-%m-%d %H:%M")}',
        '',
        'Traceback:',
        traceback_str,
        '',
        'Fix the error and restart run_pipeline.py.',
        'Completed steps will be skipped automatically.',
    ])
    send_email(EMAIL_TO, subject, body)
    print(f'[email_utils] Error email sent for step {step_num}')
