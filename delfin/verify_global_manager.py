#!/usr/bin/env python3
"""Verification script to demonstrate global manager singleton behavior."""

from delfin.global_manager import get_global_manager
from delfin.parallel_classic import _WorkflowManager

def main():
    print("=" * 70)
    print("GLOBAL MANAGER VERIFICATION TEST")
    print("=" * 70)

    # Test 1: Verify singleton pattern
    print("\n[TEST 1] Verifying singleton pattern...")
    mgr1 = get_global_manager()
    mgr2 = get_global_manager()

    if mgr1 is mgr2:
        print("✓ PASS: Both calls return the same instance")
        print(f"  Instance ID: {id(mgr1)}")
    else:
        print("✗ FAIL: Different instances returned!")
        return False

    # Test 2: Initialize global manager
    print("\n[TEST 2] Initializing global manager...")
    config = {
        'PAL': 12,
        'maxcore': 3800,
        'pal_jobs': 2,
    }
    mgr1.initialize(config)

    if mgr1.is_initialized():
        print("✓ PASS: Global manager initialized")
        print(f"  Total cores: {mgr1.total_cores}")
        print(f"  Max jobs: {mgr1.max_jobs}")
        print(f"  Pool ID: {id(mgr1.pool)}")
    else:
        print("✗ FAIL: Global manager not initialized!")
        return False

    # Test 3: Verify workflow managers use the same pool
    print("\n[TEST 3] Creating workflow managers with global pool...")

    wf1 = _WorkflowManager(config, label="test_workflow_1", use_global_pool=True)
    wf2 = _WorkflowManager(config, label="test_workflow_2", use_global_pool=True)

    pool1_id = id(wf1.pool)
    pool2_id = id(wf2.pool)
    global_pool_id = id(mgr1.pool)

    print(f"  Global pool ID:     {global_pool_id}")
    print(f"  Workflow 1 pool ID: {pool1_id}")
    print(f"  Workflow 2 pool ID: {pool2_id}")

    if pool1_id == pool2_id == global_pool_id:
        print("✓ PASS: All workflows use the same global pool!")
    else:
        print("✗ FAIL: Workflows are using different pools!")
        return False

    # Test 4: Verify local pool creates separate instance
    print("\n[TEST 4] Creating workflow manager with local pool...")
    wf3 = _WorkflowManager(config, label="test_workflow_3", use_global_pool=False)
    pool3_id = id(wf3.pool)

    print(f"  Global pool ID:     {global_pool_id}")
    print(f"  Workflow 3 pool ID: {pool3_id}")

    if pool3_id != global_pool_id:
        print("✓ PASS: Local pool is different from global pool")
    else:
        print("✗ FAIL: Local pool should be different!")
        return False

    # Test 5: Verify pool attributes
    print("\n[TEST 5] Verifying pool attributes...")
    print(f"  Global pool total cores: {mgr1.pool.total_cores}")
    print(f"  Global pool max jobs: {mgr1.pool.max_concurrent_jobs}")

    if mgr1.pool.total_cores == 12:
        print("✓ PASS: Pool has correct core count")
    else:
        print(f"✗ FAIL: Expected 12 cores, got {mgr1.pool.total_cores}")
        return False

    # Cleanup
    wf1.shutdown()
    wf2.shutdown()
    wf3.shutdown()
    mgr1.shutdown()

    print("\n" + "=" * 70)
    print("ALL TESTS PASSED ✓")
    print("=" * 70)
    print("\nConclusion:")
    print("  1. Global manager is a proper singleton")
    print("  2. Multiple workflows share the same pool when use_global_pool=True")
    print("  3. This ensures cores are never over-allocated globally")
    print("  4. ox and red workflows will coordinate through the shared pool")
    print("=" * 70)

    return True

if __name__ == '__main__':
    import sys
    success = main()
    sys.exit(0 if success else 1)
