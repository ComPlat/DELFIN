using System;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using System.Windows.Forms;

internal static class DelfinStarter {
    [DllImport("kernel32.dll", SetLastError = true)] private static extern bool AllocConsole();
    [STAThread]
    private static int Main(string[] args) {
        if (args.Length == 1 && args[0] == "--self-test") return 0;
        try {
            string folder = AppDomain.CurrentDomain.BaseDirectory;
            string script = Path.Combine(folder, "DELFIN.ps1");
            if (!File.Exists(script)) throw new FileNotFoundException("DELFIN.ps1 is missing. Reinstall DELFIN.");
            bool worker = args.Length == 3 && args[0] == "--dashboard-worker";
            Guid profile, connection;
            if (args.Length != 0 && !worker) throw new ArgumentException("Invalid starter arguments.");
            string extra = "";
            if (worker) {
                if (!Guid.TryParse(args[1], out profile) || !Guid.TryParse(args[2], out connection))
                    throw new ArgumentException("Invalid connection identifiers.");
                if (!AllocConsole()) throw new InvalidOperationException("Cannot open the SSH login console.");
                extra = " -Mode Dashboard -ProfileId " + profile.ToString() + " -ConnectionId " + connection.ToString();
            }
            ProcessStartInfo info = new ProcessStartInfo();
            info.FileName = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.Windows),
                                        @"System32\WindowsPowerShell\v1.0\powershell.exe");
            info.Arguments = "-NoProfile -ExecutionPolicy RemoteSigned -STA -File \"" + script + "\"" + extra;
            info.WorkingDirectory = folder;
            info.UseShellExecute = false;
            info.CreateNoWindow = !worker;
            info.RedirectStandardError = !worker;
            info.RedirectStandardOutput = !worker;
            using (Process child = Process.Start(info)) {
                if (worker) {
                    child.WaitForExit();
                    return child.ExitCode;
                }
                var errors = child.StandardError.ReadToEndAsync();
                var output = child.StandardOutput.ReadToEndAsync();
                child.WaitForExit();
                if (child.ExitCode != 0)
                    MessageBox.Show(errors.Result, "DELFIN startup failed", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return child.ExitCode;
            }
        } catch (Exception error) {
            MessageBox.Show(error.Message, "DELFIN", MessageBoxButtons.OK, MessageBoxIcon.Error);
            return 1;
        }
    }
}
