using System;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using System.Threading;
using System.Security.Principal;

internal static class DelfinStarter {
    [DllImport("kernel32.dll", SetLastError = true)] private static extern bool AllocConsole();
    [DllImport("shell32.dll", CharSet=CharSet.Unicode)] private static extern int SetCurrentProcessExplicitAppUserModelID(string id);
    [DllImport("user32.dll", CharSet=CharSet.Unicode)] private static extern IntPtr FindWindow(string className, string title);
    [DllImport("user32.dll")] private static extern bool ShowWindow(IntPtr window, int command);
    [DllImport("user32.dll")] private static extern bool SetForegroundWindow(IntPtr window);
    private static Mutex appMutex;
    private static bool ownsMutex;
    [STAThread]
    private static int Main(string[] args) {
        if (args.Length == 1 && args[0] == "--self-test") return 0;
        try {
            string folder = AppDomain.CurrentDomain.BaseDirectory;
            string script = Path.Combine(folder, "DELFIN.ps1");
            if (!File.Exists(script)) throw new FileNotFoundException("DELFIN.ps1 is missing. Reinstall DELFIN.");
            bool smoke = args.Length == 1 && args[0] == "--gui-smoke-test";
            bool worker = args.Length == 3 && args[0] == "--dashboard-worker";
            if (!worker && !smoke) {
                bool created;
                string sid = WindowsIdentity.GetCurrent().User.Value;
                appMutex = new Mutex(true, @"Local\DELFIN_Launcher_" + sid, out created);
                ownsMutex = created;
                if (!created) {
                    IntPtr window = FindWindow(null, "DELFIN - SSH Dashboard");
                    if (window != IntPtr.Zero) { ShowWindow(window, 9); SetForegroundWindow(window); }
                    else MessageBox.Show("DELFIN is already running or starting. Use its existing settings window.", "DELFIN");
                    return 0;
                }
            }
            Guid profile, connection;
            if (args.Length != 0 && !worker && !smoke) throw new ArgumentException("Invalid starter arguments.");
            string extra = smoke ? " -GuiSmokeTest" : "";
            if (worker) {
                if (!Guid.TryParse(args[1], out profile) || !Guid.TryParse(args[2], out connection))
                    throw new ArgumentException("Invalid connection identifiers.");
                SetCurrentProcessExplicitAppUserModelID("ComPlat.DELFIN.SSH");
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
        } finally {
            if (appMutex != null) {
                if (ownsMutex) appMutex.ReleaseMutex();
                appMutex.Dispose();
            }
        }
    }
}
