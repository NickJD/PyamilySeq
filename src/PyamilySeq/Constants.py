import subprocess

PyamilySeq_Version = 'v0.3.0'



def is_tool_installed(tool_name):
    """Check if a tool is installed and available in PATH."""
    try:
        subprocess.run([tool_name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except subprocess.CalledProcessError:
        return False
    except FileNotFoundError:
        return False