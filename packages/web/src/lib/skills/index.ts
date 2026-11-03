import contributeContent from '../../../content/skills/omnipath-contribute/SKILL.md?raw';
import parquetContent from '../../../content/skills/omnipath-parquet/SKILL.md?raw';

export type Skill = {
  name: string;
  title: string;
  /** What the skill is for, using the data or contributing to it; it picks the icon. */
  kind: 'Use the data' | 'Contribute';
  description: string;
  href: string;
  content: string;
};

export const skills: Skill[] = [
  {
    name: 'omnipath-parquet',
    title: 'Use OmniPath',
    description:
      'Install resources with the OmniPath client and query entities and relations locally.',
    kind: 'Use the data',
    href: '/skills/omnipath-parquet/SKILL.md',
    content: parquetContent,
  },
  {
    name: 'omnipath-contribute',
    title: 'Contribute to OmniPath',
    description: 'Teach your agent to add a data resource and submit it to pypath.',
    kind: 'Contribute',
    href: '/skills/omnipath-contribute/SKILL.md',
    content: contributeContent,
  },
];

export function getSkill(name: string) {
  return skills.find((skill) => skill.name === name);
}
