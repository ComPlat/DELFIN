(() => {
  'use strict';
  const page = document.querySelector('.framework-page');
  if (!page) return;
  const components = {
    person: ['GOAL & SCOPE', 'The person', 'Defines the task, working scope and decisions. Approvals come from the person through the applicable approval route; a model-generated plan does not authorise itself.'],
    context: ['RETRIEVED INPUT', 'Context & memory', 'Project instructions, relevant facts, invoked skills and retrieved history help assemble the next request. Retrieval selects information within a budget; a memory is background data, not permission.'],
    engine: ['CONVERSATION LOOP', 'Engine & model', 'The engine assembles the conversation and task state; the backend connects to the selected model. A proposed tool call crosses the permission boundary before execution. Observed results then feed the next round.'],
    state: ['CONTINUITY', 'Task & session state', 'Active goals, messages and working state keep the session focused. Compaction bounds the context while retaining recent rounds and important outcomes; durable memories and indexed history have separate lifetimes.'],
    gate: ['ALLOW / ASK / REFUSE', 'Permission boundary', 'The proposed action is evaluated before execution. Path rules, command rules, protected routes and host-dependent isolation determine what may run, what needs approval and what is refused.'],
    tools: ['EXECUTION', 'Tools & subagents', 'File, code, document and chemistry tools perform the work. Subagents carry separate contexts and budgets with inherited boundaries. A requested action and its actual execution are distinct events.'],
    evidence: ['OBSERVED RESULTS', 'Evidence & checks', 'Tool results, exit status, tests and calculation files supply evidence. Scientific output diagnostics and supported claim checks ask different questions; coverage depends on the route and available records.'],
    answer: ['GROUNDED HANDOVER', 'Answer & next step', 'The response reports the deliverable, observed results and verification limits. Supported completion and claim checks can flag assertions that outrun the recorded work. The next action remains subject to the same boundary.']
  };
  function describe(panel, eyebrow, title, text) {
    panel.querySelector('.eyebrow').textContent = eyebrow;
    panel.querySelector('h3').textContent = title;
    panel.querySelector('p').textContent = text;
  }
  const componentButtons = [...page.querySelectorAll('[data-component]')];
  const componentPanel = page.querySelector('#system-detail');
  if (componentPanel) componentButtons.forEach((button) => {
    const component = components[button.dataset.component];
    if (!component) return;
    button.disabled = false;
    button.addEventListener('click', () => {
      componentButtons.forEach((other) => other.setAttribute('aria-pressed', String(other === button)));
      describe(componentPanel, ...component);
    });
  });
  const taskButtons = [...page.querySelectorAll('[data-task-step]')];
  const taskSteps = [...page.querySelectorAll('#task-step-source > li')];
  const taskPanel = page.querySelector('#task-step-detail');
  function selectStep(button) {
    const index = Number(button.dataset.taskStep);
    const source = taskSteps[index];
    if (!source || !taskPanel) return;
    taskButtons.forEach((other) => other.setAttribute('aria-pressed', String(other === button)));
    describe(taskPanel, String(index + 1).padStart(2, '0') + ' / ' + button.textContent.slice(2).trim().toUpperCase(),
      source.querySelector('h3').textContent, source.querySelector('p').textContent);
  }
  taskButtons.forEach((button) => {
    if (!taskSteps[Number(button.dataset.taskStep)] || !taskPanel) return;
    button.disabled = false;
    button.addEventListener('click', () => selectStep(button));
  });
  if (taskButtons[0]) selectStep(taskButtons[0]);
  function openCheck(hash) {
    const check = [...page.querySelectorAll('.framework-checks details')].find((node) => '#' + node.id === hash);
    if (check) check.open = true;
  }
  page.querySelectorAll('.validation-track a').forEach((link) => {
    link.addEventListener('click', () => openCheck(link.hash));
  });
  addEventListener('hashchange', () => openCheck(location.hash));
  openCheck(location.hash);
})();
