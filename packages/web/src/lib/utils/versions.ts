/** Ascending order, matching the API's numeric-version / legacy-time selection. */
export function compareResourceVersions(a: string, b: string, aTime = 0, bTime = 0): number {
  const numeric = /^[0-9]+(?:\.[0-9]+)*$/;
  const aNumeric = numeric.test(a);
  const bNumeric = numeric.test(b);
  if (aNumeric !== bNumeric) return aNumeric ? 1 : -1;
  if (aNumeric) {
    const left = a.split('.').map(BigInt);
    const right = b.split('.').map(BigInt);
    for (let i = 0; i < Math.min(left.length, right.length); i++) {
      if (left[i] !== right[i]) return left[i] > right[i] ? 1 : -1;
    }
    if (left.length !== right.length) return left.length - right.length;
  } else if (aTime !== bTime) {
    return aTime - bTime;
  }
  return a === b ? 0 : a > b ? 1 : -1;
}
