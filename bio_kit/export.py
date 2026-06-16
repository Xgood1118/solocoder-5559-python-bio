from .sequence_io import write_fasta, write_genbank

def export_fasta(sequences, file_path):
    write_fasta(sequences, file_path)
    return file_path

def export_genbank(sequences, file_path):
    write_genbank(sequences, file_path)
    return file_path

def export_svg(svg_content, file_path):
    with open(file_path, 'w') as f:
        f.write(svg_content)
    return file_path

def export_newick(tree, file_path):
    newick = tree.to_newick()
    with open(file_path, 'w') as f:
        f.write(newick)
    return file_path

def export_png(svg_content, file_path):
    try:
        from cairosvg import svg2png
        svg2png(bytestring=svg_content.encode(), write_to=file_path)
        return file_path
    except ImportError:
        try:
            import base64
            html_content = f'''
            <!DOCTYPE html>
            <html>
            <head><title>Export</title></head>
            <body>{svg_content}</body>
            </html>
            '''
            alt_path = file_path.replace('.png', '.svg')
            with open(alt_path, 'w') as f:
                f.write(svg_content)
            return alt_path
        except:
            pass
    
    alt_path = file_path.replace('.png', '.svg')
    with open(alt_path, 'w') as f:
        f.write(svg_content)
    return alt_path
