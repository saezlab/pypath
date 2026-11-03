/** How a relation reads between its two entities: "Increases activity of", "via binding". */
export type RelationPhrase = {
  text: string;
  /** The phrase for a narrow layout: the arrow says the direction, so "phosphorylation". */
  short: string;
  /** How the effect happens, e.g. "via binding". */
  detail?: string;
  direction?: 'up' | 'down';
};

export type PhraseRelation = {
  predicate?: string | null;
  sign?: number | null;
  displayLabel?: string | null;
  qualifiers?: Record<string, string> | null;
};

const humanize = (value: string) => value.replace(/_/g, ' ').toLowerCase();

// Mechanisms already said by a qualified label ("Binds", "Transports").
const IMPLIED_MECHANISM: Record<string, string> = {
  Binds: 'binding',
  Transports: 'relocalization',
};

/** `fallbackLabel` is the plain predicate label, used when there is nothing more specific. */
export function relationPhrase(relation: PhraseRelation, fallbackLabel: string): RelationPhrase {
  const qualifiers = relation.qualifiers ?? {};
  const predicate = (relation.predicate ?? '').replace(/^biolink:/i, '').toLowerCase();
  const directionValue =
    qualifiers.object_direction_qualifier ??
    (predicate === 'affects' && relation.sign === 1
      ? 'increased'
      : predicate === 'affects' && relation.sign === -1
        ? 'decreased'
        : undefined);
  const direction =
    directionValue === 'increased' ? 'up' : directionValue === 'decreased' ? 'down' : undefined;
  const mechanism = qualifiers.causal_mechanism_qualifier;
  const label = relation.displayLabel?.trim();
  const detail =
    mechanism && (!label || IMPLIED_MECHANISM[label] !== mechanism)
      ? `via ${humanize(mechanism)}`
      : undefined;

  if (label) return { text: label, short: label, detail, direction };
  if (predicate === 'affects') {
    const verb = direction === 'up' ? 'Increases' : direction === 'down' ? 'Decreases' : 'Affects';
    const aspect = qualifiers.object_aspect_qualifier;
    if (!aspect) return { text: verb, short: verb, detail, direction };
    return {
      text: `${verb} ${humanize(aspect)} of`,
      short: direction ? humanize(aspect) : `Affects ${humanize(aspect)}`,
      detail,
      direction,
    };
  }
  return { text: fallbackLabel, short: fallbackLabel, detail, direction };
}
