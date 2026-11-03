import { getSkill } from '$lib/skills';
import { error, type RequestHandler } from '@sveltejs/kit';

export const GET: RequestHandler = ({ params }) => {
  const skill = getSkill(params.skill ?? '');
  if (!skill) error(404, 'Skill not found');

  return new Response(skill.content, {
    headers: {
      'Content-Type': 'text/markdown; charset=utf-8',
      'Content-Disposition': `attachment; filename="${skill.name}-SKILL.md"`,
    },
  });
};
