"""Native layout, selection and scroll regression checks for the chat renderer."""
import pytest

from delfin.dashboard import tab_agent as T
from test_the_chat_followed_a_reader_who_had_scrolled_away import _chat_scroll_js

pytestmark = pytest.mark.slow


@pytest.fixture
def page():
    api = pytest.importorskip('playwright.sync_api')
    with api.sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={'width': 1000, 'height': 700})
        page.set_content(T._AGENT_CSS + '''
          <div class="delfin-agent-chat-host" style="height:400px;width:800px">
            <div class="widget-html-content"></div>
          </div>''')
        page.add_script_tag(content=_chat_scroll_js())
        yield page
        browser.close()


def _show(page, content, *, streaming=True):
    rendered = T._streaming_html(content) if streaming else T._md_to_html(content)
    value = T._chat_value('<div class="delfin-agent-chat">'
                          '<div class="delfin-chat-msg delfin-chat-agent">'
                          '<div class="delfin-chat-role">Agent</div>' + rendered +
                          '</div></div>')
    page.evaluate('value => {document.querySelector(".widget-html-content").innerHTML=value}', value)
    page.wait_for_function('''() => document.querySelector('.widget-html-content img') === null''')


def test_empty_chat_has_no_phantom_scroll(page):
    _show(page, 'Start a conversation...', streaming=False)
    assert page.locator('.delfin-agent-chat-host').evaluate(
        '(c) => c.scrollTop === 0 && c.scrollHeight === c.clientHeight')


def test_finished_paragraph_selection_survives_tokens_and_finalization(page):
    head = '**Read this paragraph** while more text arrives.\n\n'
    _show(page, head + 'A')
    page.evaluate('''() => {
      window.keptParagraph = document.querySelector('.delfin-chat-agent p');
      const range = document.createRange(); range.selectNodeContents(keptParagraph);
      getSelection().removeAllRanges(); getSelection().addRange(range);
      window.keptSelection = getSelection().toString();
    }''')
    for tail in ['A few', 'A few more tokens.', 'A few more tokens.\n\nNext paragraph']:
        _show(page, head + tail)
        assert page.evaluate('''() => keptParagraph === document.querySelector('.delfin-chat-agent p')
            && getSelection().toString() === keptSelection''')
    _show(page, head + 'A few more tokens.\n\nNext paragraph', streaming=False)
    assert page.evaluate('() => getSelection().toString() === keptSelection')


def test_long_reply_does_not_reflow_or_move_a_reader_at_the_old_limit(page):
    paragraph = '## Result\n\n**Stable** paragraph with `code`.\n\n'
    before = paragraph * (19900 // len(paragraph))
    _show(page, before + 'tail')
    chat = page.locator('.delfin-agent-chat-host')
    chat.hover()
    page.mouse.wheel(0, -100000)
    page.wait_for_function('() => document.querySelector(".delfin-agent-chat-host").scrollTop === 0')
    page.evaluate('''() => {
      window.keptHeading = document.querySelector('.delfin-chat-agent h2');
      window.headingY = keptHeading.getBoundingClientRect().top;
      window.initialHeight = document.querySelector('.delfin-agent-chat-host').scrollHeight;
    }''')
    _show(page, before + paragraph * 10 + 'tail')
    assert page.evaluate('''() => {
      const c = document.querySelector('.delfin-agent-chat-host');
      return c.scrollTop === 0 && keptHeading === document.querySelector('.delfin-chat-agent h2')
        && keptHeading.getBoundingClientRect().top === headingY
        && c.scrollHeight - initialHeight < 1000 && c.dataset.delfinFollow === '0';
    }''')
