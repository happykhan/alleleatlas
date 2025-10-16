"""CLI wrapper for the alleleatlas pipeline.

This is a small, user-friendly CLI that wraps `alleleatlas.main.main` and
prints a neat summary of the invocation before running.
"""

from pathlib import Path
import typer
from rich.console import Console
from rich.table import Table

from alleleatlas.main import main

app = typer.Typer(help='Run the alleleatlas cgMLST -> HierCC -> embeddings pipeline')
console = Console()


def _print_invocation(profile: Path, outdir: Path, nproc: int, force: bool, skip_cluster: bool, skip_embed: bool):
    tbl = Table(title='alleleatlas invocation')
    tbl.add_column('option', style='cyan')
    tbl.add_column('value', style='magenta')
    tbl.add_row('profile', str(profile))
    tbl.add_row('outdir', str(outdir))
    tbl.add_row('nproc', str(nproc))
    tbl.add_row('force', str(force))
    tbl.add_row('skip_cluster', str(skip_cluster))
    tbl.add_row('skip_embed', str(skip_embed))
    console.print(tbl)


@app.command()
def run(
    profile: Path = typer.Argument(..., help='Path to a cgMLST profile file (tsv/csv/.gz/.xz)'),
    outdir: Path = typer.Argument('output', help='Output directory (default: output)'),
    nproc: int = typer.Option(2, '-n', '--nproc', help='Number of worker processes'),
    force: bool = typer.Option(False, '-f', '--force', help='Recompute all steps, ignore cache'),
    skip_cluster: bool = typer.Option(False, '--skip-cluster', help='Skip HierCC clustering and evaluation'),
    skip_embed: bool = typer.Option(False, '--skip-embed', help='Skip UMAP/MDS embedding step'),
):
    """Run the full alleleatlas pipeline with caching support.

    Uses cache to skip expensive steps if outputs exist (unless --force is given).
    """
    # Ensure outdir exists
    outdir.mkdir(parents=True, exist_ok=True)

    _print_invocation(profile, outdir, nproc, force, skip_cluster, skip_embed)

    if skip_cluster or skip_embed or nproc != 2:
        console.print('[yellow]Note: skip_cluster and skip_embed flags are not yet wired to internal steps.[/yellow]')

    console.print(f'[green]Starting pipeline[/green] (force={force}): {profile} -> {outdir}')
    try:
        main(str(profile), str(outdir), force=force, nproc=nproc)
        console.print('[bold green]✓ Pipeline finished successfully.[/bold green]')
    except Exception as e:
        console.print(f'[bold red]✗ Pipeline failed:[/bold red] {e}')
        raise


if __name__ == '__main__':
    app()
