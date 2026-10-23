using System;
using System.Runtime.InteropServices;

public static class DelfinTaskbar {
    [StructLayout(LayoutKind.Sequential)] private struct PropertyKey {
        public Guid Format;
        public uint Id;
        public PropertyKey(uint id) { Format = new Guid("9F4C2855-9F79-4B39-A8D0-E1D42DE1D5F3"); Id = id; }
    }
    // PROPVARIANT is 24 bytes on x64 and 16 bytes on x86.
    [StructLayout(LayoutKind.Sequential)] private struct Variant {
        public ushort Type, Reserved1, Reserved2, Reserved3;
        public IntPtr Pointer, Padding;
    }
    [ComImport, Guid("886D8EEB-8CF2-4446-8D02-CDBA1DBDCF99"), InterfaceType(ComInterfaceType.InterfaceIsIUnknown)]
    private interface IPropertyStore {
        [PreserveSig] int GetCount(out uint count);
        [PreserveSig] int GetAt(uint index, out PropertyKey key);
        [PreserveSig] int GetValue(ref PropertyKey key, out Variant value);
        [PreserveSig] int SetValue(ref PropertyKey key, ref Variant value);
        [PreserveSig] int Commit();
    }
    [DllImport("shell32.dll")] private static extern int SHGetPropertyStoreForWindow(IntPtr window, ref Guid iid, [MarshalAs(UnmanagedType.Interface)] out IPropertyStore store);
    [DllImport("ole32.dll")] private static extern int PropVariantClear(ref Variant value);
    private static IPropertyStore Store(IntPtr window) {
        Guid iid = typeof(IPropertyStore).GUID;
        IPropertyStore store;
        Marshal.ThrowExceptionForHR(SHGetPropertyStoreForWindow(window, ref iid, out store));
        return store;
    }
    private static void Set(IPropertyStore store, uint id, string text) {
        PropertyKey key = new PropertyKey(id);
        Variant value = new Variant();
        value.Type = 31; // VT_LPWSTR
        value.Pointer = Marshal.StringToCoTaskMemUni(text);
        try { Marshal.ThrowExceptionForHR(store.SetValue(ref key, ref value)); }
        finally { PropVariantClear(ref value); }
    }
    public static void Configure(IntPtr window, string starter, string icon) {
        IPropertyStore store = Store(window);
        try {
            Set(store, 2, "\"" + starter + "\"");
            Set(store, 3, icon + ",0");
            Set(store, 4, "DELFIN");
            Set(store, 5, "ComPlat.DELFIN.Launcher");
        } finally { Marshal.ReleaseComObject(store); }
    }
    public static string Read(IntPtr window, uint id) {
        IPropertyStore store = Store(window);
        Variant value = new Variant();
        try {
            PropertyKey key = new PropertyKey(id);
            Marshal.ThrowExceptionForHR(store.GetValue(ref key, out value));
            return value.Type == 31 ? Marshal.PtrToStringUni(value.Pointer) : "";
        } finally { PropVariantClear(ref value); Marshal.ReleaseComObject(store); }
    }
}
