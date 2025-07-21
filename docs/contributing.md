# Contributing

## Development Setup

### Environment Setup

```bash
# Clone the repository
git clone https://github.com/poclab-web/streamlit-pyscf.git
cd streamlit-pyscf

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install development dependencies
pip install -r docs/requirements.txt
```

### Running Tests

```bash
# Run unit tests
python -m pytest tests/

# Run with coverage
python -m pytest --cov=logic tests/
```

### Building Documentation

```bash
# Build HTML documentation
cd docs
sphinx-build -b html . _build/html

# Build and serve locally
sphinx-autobuild . _build/html
```

## Contribution Guidelines

### Code Style

- Follow PEP 8 style guide
- Use type hints where possible
- Add docstrings to all functions and classes
- Maximum line length: 88 characters

### Git Workflow

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new features
5. Update documentation
6. Submit a pull request

### Pull Request Process

1. Ensure all tests pass
2. Update documentation if needed
3. Add a clear description of changes
4. Reference any related issues

## Code Structure

### Adding New Calculations

1. Add core logic to `logic/` directory
2. Create corresponding page in `pages/`
3. Update configuration if needed
4. Add tests and documentation

### Database Schema Changes

When modifying the database schema:

1. Update `logic/database.py`
2. Provide migration script if needed
3. Update documentation

### Adding New Dependencies

1. Add to `requirements.txt`
2. Update `docs/requirements.txt` if needed
3. Test compatibility

## Issue Reporting

When reporting issues, please include:

- Python version
- Operating system
- Error messages (full traceback)
- Minimal example to reproduce
- Expected vs actual behavior

## Development Priorities

### High Priority
- Performance optimization
- Error handling improvements
- User interface enhancements

### Medium Priority
- Additional calculation methods
- Export/import features
- Batch processing capabilities

### Future Plans
- Web API development
- Cloud deployment support
- Integration with other quantum chemistry packages
