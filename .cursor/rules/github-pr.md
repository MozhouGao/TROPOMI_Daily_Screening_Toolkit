# GitHub Pull Request Rules

## PR Creation Process

1. **Always use the template**: Use `.github/pull_request_template.md` as the base template
2. **Check current branch and git status** before creating PR
3. **Review recent commits** to understand changes
4. **Generate descriptive title** following conventional commit format
5. **Fill out all template sections** completely

## PR Title Format

- Follow conventional commit format: `<type>: <description>`
- **Type**: Use sentence case (only first word capitalized)
- **Description**: Use lowercase (except for proper nouns, acronyms, or brand names)
- Keep under 100 characters
- Examples:
  - "Feat: add user authentication with Auth0"
  - "Refactor: optimize asset loading performance"
  - "Fix: resolve map rendering issues"

## PR Body Structure

Use the template structure from `.github/pull_request_template.md`:

### Required Sections:
- **✍️ Description**: Clear summary of changes and motivation
- **📝 Type of Change**: Select appropriate checkboxes
- **🔍 What Changed**: Bullet points of key changes
- **📷 Screenshots**: If UI changes are involved
- **🔗 Related Issues**: Link to issues, discussions, or documentation
- **📝 Additional Notes**: Any extra context for reviewers

### Optional Sections:
- **🧪 Testing**: Testing checklist (if not already covered in template)
- **📋 Checklist**: Quality assurance items

## Template Usage Guidelines

- **Fill out all sections** - don't leave empty sections
- **Be specific** about what changed and why
- **Include screenshots** for UI changes
- **Link related issues** using GitHub keywords (Closes #, Fixes #, etc.)
- **Check all applicable boxes** in Type of Change section
- **Provide clear testing information** if manual testing was performed

## Example PR Body Structure

```markdown
## ✍️ Description

Enhanced Loading component with logo animation and improved styling. Replaced Ant Design Spin component with custom Loading component for consistent loading experience across the application.

## 📝 Type of Change

- [x] 🚀 New feature (adds functionality)
- [x] 🧹 Refactor (no functional changes)

## 🔍 What Changed

- Enhanced Loading component with logo animation
- Replaced Ant Design Spin with custom Loading component
- Improved visual feedback during loading states
- Added consistent styling across loading scenarios

## 📷 Screenshots (if applicable)

[Add screenshots here if UI changes were made]

## 🔗 Related Issues

- Closes #123
- Related to #456

## 📝 Additional Notes

No breaking changes introduced. All existing functionality preserved.
```

**Example PR Title:** `feat: enhance loading component with logo animation`

## Quality Checklist

Before creating PR:
- [ ] All changes are committed and pushed
- [ ] Code follows project style guidelines
- [ ] Tests pass (if applicable)
- [ ] No console.log statements in production code
- [ ] Internationalization keys added (if applicable)
- [ ] Custom Loading component used instead of Ant Design Spin
- [ ] Self-review completed
